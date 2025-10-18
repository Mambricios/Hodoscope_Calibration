#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <algorithm>
#include <limits>
#include <numeric>

// Declaraciones de funciones
bool safe_stof(const std::string& str, float& result);
void CrearAchivoROOT();
void GenerarHistograma(std::vector<std::vector<float>>& all_signal_sums);
std::vector<float> gaussianas();
void num_phe_candidatos(const std::vector<float>& mu_values, const std::vector<std::vector<float>>& all_signal_sums);
void ProcesarAmplitudesMaximas(const std::vector<float>& mu_values);
void ProcesarDatosFotoelectrones(const std::vector<float>& mu_values, std::vector<std::vector<float>>& all_signal_sums);

bool safe_stof(const std::string& str, float& result) {
try {
result = std::stof(str);
return true;
} catch (const std::invalid_argument& e) {
std::cerr << "Error: valor no numérico detectado: " << str << std::endl;
return false;
} catch (const std::out_of_range& e) {
std::cerr << "Error: valor fuera de rango detectado: " << str << std::endl;
return false;
}
}

void CrearAchivoROOT() {
// Constantes
const float bin_width = 8.0;
const int bins_per_record = 1024;
const int num_channels = 32;
const int min_channel = 0;
const int max_channel = 31;

// Abrir archivo CSV

std::ifstream archivo("500LinesExample.csv");

if (!archivo.is_open()) {

    std::cerr << "No se pudo abrir el archivo." << std::endl;

    return;

}


// Crear archivo ROOT

TFile *file = new TFile("proyecto2_funciones2.root", "RECREATE");


// Crear un vector de TTree para cada canal

std::vector<TTree*> trees(num_channels);

std::vector<float> Data(bins_per_record);

float Time;

int Bin;

int Bool;

float Offset;


for (int i = 0; i < num_channels; ++i) {

    trees[i] = new TTree(Form("Channel%d", i), Form("Datos del canal %d", i));

    trees[i]->Branch("Time", &Time, "Time/F");

    trees[i]->Branch("Bin", &Bin, "Bin/I");

    trees[i]->Branch("Bool", &Bool, "Bool/I");

    trees[i]->Branch("Offset", &Offset, "Offset/F");

    trees[i]->Branch("Data", &Data[0], "Data[1024]/F");

}


std::string linea;


while (std::getline(archivo, linea)) {

    std::vector<std::string> columnas;

    std::stringstream ss(linea);

    std::string campo;


    while (std::getline(ss, campo, ';')) {

        columnas.push_back(campo);

    }


    if (columnas.size() >= 1029) {

        if (!safe_stof(columnas[0], Time)) continue;

        if (!safe_stof(columnas[4], Offset)) continue;

        int Channel = std::stoi(columnas[1]);


        if (Channel < min_channel || Channel > max_channel) {

            std::cerr << "Canal " << Channel << " fuera de rango, ignorando fila." << std::endl;

            continue;

        }


        Bin = std::stoi(columnas[2]);

        Bool = std::stoi(columnas[3]);


        for (int i = 0; i < bins_per_record; ++i) {

            if (!safe_stof(columnas[5 + i], Data[i])) continue;

        }


        trees[Channel]->Fill();

    }

}


// Escribir los árboles en el archivo ROOT

for (int i = 0; i < num_channels; ++i) {

    trees[i]->Write();

    delete trees[i];

}


file->Close();

delete file;

}

void GraficarDatosReajustados(std::vector<std::vector<float>>& all_signal_sums) {
    const int num_channels = 32;
    const float bin_width = 8.0;
    const int bins_per_record = 1024;

    TFile *file = new TFile("proyecto2_funciones2.root", "UPDATE");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: No se pudo abrir el archivo ROOT." << std::endl;
        return;
    }

    for (int channel = 0; channel < num_channels; ++channel) {
        TTree* tree = (TTree*)file->Get(Form("Channel%d", channel));
        float Data[bins_per_record];
        tree->SetBranchAddress("Data", Data);

        // Histograma 1D original
        TH1F* hist_reajustado = new TH1F(Form("hist_reajustado_Channel%d", channel),
                                        Form("Adjusted Data Channel %d", channel),
                                        bins_per_record, 0, bins_per_record * bin_width);
        hist_reajustado->GetXaxis()->SetTitle("Time (ns)");
        hist_reajustado->GetYaxis()->SetTitle("Amplitude, ADC");

        // Histograma de persistencia (2D)
        TH2F* hist_persistencia = new TH2F(Form("hist_persistence_Channel%d", channel),
                                           Form("Persistence Histogram Channel %d", channel),
                                           bins_per_record, 0, bins_per_record * bin_width,
                                           400, -200, 800);  // rango Y ajustable
        hist_persistencia->GetXaxis()->SetTitle("Time (ns)");
        hist_persistencia->GetYaxis()->SetTitle("Amplitude (ADC)");
        hist_persistencia->GetZaxis()->SetTitle("Frecuencia (Persistencia)");

        int nEntries = tree->GetEntries();

        for (int entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);

            // Calcular baseline
            float baseline = 0.0;
            int n_baseline = 300;
            for (int i = 0; i < n_baseline; ++i) baseline += Data[i];
            baseline /= n_baseline;

            // Reajustar datos (invertir y restar baseline)
            for (int i = 0; i < bins_per_record; ++i) {
                Data[i] = -(Data[i] - baseline);
            }

            // Llenar histograma 1D (promedio del último evento)
            for (int i = 0; i < bins_per_record; ++i) {
                hist_reajustado->SetBinContent(i + 1, Data[i]);
            }

            // Llenar histograma de persistencia
            for (int i = 0; i < bins_per_record; ++i) {
                float time_ns = i * bin_width;
                hist_persistencia->Fill(time_ns, Data[i]);
            }
        }

        // Guardar ambos histogramas
        hist_reajustado->Write();
        hist_persistencia->Write();

        delete hist_reajustado;
        delete hist_persistencia;
    }

    file->Close();
    delete file;
}


void GenerarHistograma(std::vector<std::vector<float>>& all_signal_sums) {
const int num_channels = 32;
const int bins_per_record = 1024;
const int bin_width = 8.0;

TFile *file = new TFile("proyecto2_funciones2.root", "UPDATE");

TH1F* hist_signal_sum[num_channels];


for (int i = 0; i < num_channels; ++i) {

    hist_signal_sum[i] = new TH1F(Form("hist_signal_sum_Channel%d", i),

                                   Form("Charge Distribution Channel %d", i),

                                   100, 0, 3000);

    hist_signal_sum[i]->GetXaxis()->SetTitle("ADC*ns");

    hist_signal_sum[i]->GetYaxis()->SetTitle("Entries");

}


for (int channel = 0; channel < num_channels; ++channel) {

    TTree* tree = (TTree*)file->Get(Form("Channel%d", channel));

    float Data[bins_per_record];

    tree->SetBranchAddress("Data", Data);


    int nEntries = tree->GetEntries();


    for (int entry = 0; entry < nEntries; ++entry) {

        tree->GetEntry(entry);


        float baseline = 0.0;

        int n_baseline = 300;

        for (int i = 0; i < n_baseline; ++i) baseline += Data[i];

        baseline /= n_baseline;


        TH1F *hist_temp = new TH1F("hist_temp", "Histograma Temporal para Baseline", n_baseline, 0, n_baseline);

        for (int i = 0; i < n_baseline; ++i) {

            hist_temp->Fill(Data[i] - baseline);

        }


        float sigma = hist_temp->GetRMS();

        delete hist_temp;


        for (int i = 0; i < bins_per_record; ++i) {

            Data[i] = -(Data[i] - baseline);

        }


        std::vector<float> signal_sum;

        for (int i = 1; i < bins_per_record - 5; ++i) {

            float signal_amplitude = Data[i];

            if (signal_amplitude >= 5 * sigma && signal_amplitude <= 20 * sigma) {

                float sum = 0.0;

                for (int j = i - 2; j <= i + 4; ++j) {

                    sum += Data[j];

                }
                sum *= bin_width;

                signal_sum.push_back(sum);

                i += 4;

            }

        }


        // Almacenar las sumas de señales

        for (float sum : signal_sum) {

            all_signal_sums[channel].push_back(sum);

            hist_signal_sum[channel]->Fill(sum);

        }

    }

}


// Escribir los histogramas

for (int i = 0; i < num_channels; ++i) {

    hist_signal_sum[i]->Write();

    delete hist_signal_sum[i];

}


file->Close();

delete file;

}

std::vector<float> gaussianas() {
const int num_channels = 32;
std::vector<float> mu_values(num_channels, 0.0);
TSpectrum spectrum;

TFile *file = new TFile("proyecto2_funciones2.root", "UPDATE");

TH1F* hist_signal_sum[num_channels];

TCanvas* canvases[num_channels];


for (int i = 0; i < num_channels; ++i) {

    hist_signal_sum[i] = (TH1F*)file->Get(Form("hist_signal_sum_Channel%d", i));

    canvases[i] = new TCanvas(Form("canvas_ch%d", i), Form("Channel %d Peak Analysis", i), 800, 600);

}


for (int i = 0; i < num_channels; ++i) {

    canvases[i]->cd();

    hist_signal_sum[i]->Draw();


    int n_peaks = spectrum.Search(hist_signal_sum[i], 1, "nobackground", 0.1);

    hist_signal_sum[i]->Draw();


    if (n_peaks > 0) {

        Double_t *xpeaks = spectrum.GetPositionX();

        float closest_peak = 0.0;

        float min_distance = std::numeric_limits<float>::max();


        for (int j = 0; j < n_peaks; ++j) {

            float distance = std::abs(xpeaks[j] - 700.0);

            if (distance < min_distance) {

                min_distance = distance;

                closest_peak = xpeaks[j];

            }

        }


        float fit_range = 50.0;

        float min_fit = std::max<double>(closest_peak - fit_range / 2.0, hist_signal_sum[i]->GetXaxis()->GetXmin());

        float max_fit = std::min<double>(closest_peak + fit_range / 2.0, hist_signal_sum[i]->GetXaxis()->GetXmax());

        float height = hist_signal_sum[i]->GetBinContent(hist_signal_sum[i]->FindBin(closest_peak));

        float sigma_est = hist_signal_sum[i]->GetRMS();


        TF1 *gaus_fit = new TF1(Form("gaus_spe_%d", i), "gaus", min_fit, max_fit);

        gaus_fit->SetParameters(height, closest_peak, sigma_est);

        gaus_fit->SetLineColor(kRed);

        hist_signal_sum[i]->Fit(gaus_fit, "RQ+");

        mu_values[i] = gaus_fit->GetParameter(1);


        TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

        leg->AddEntry(hist_signal_sum[i], "Datos", "l");

        leg->AddEntry(gaus_fit, "Ajuste Gaussiano", "l");

        leg->Draw();


        canvases[i]->Write(Form("peak_analysis_ch%d", i), TObject::kOverwrite);


        delete gaus_fit;

        delete leg;

    } else {

        std::cerr << "Advertencia: No se encontraron picos en el canal " << i << std::endl;

        hist_signal_sum[i]->SetTitle(Form("Channel %d - No peaks found", i));

        canvases[i]->Write(Form("peak_analysis_ch%d", i), TObject::kOverwrite);

    }

}


for (int i = 0; i < num_channels; ++i) {

    delete hist_signal_sum[i];

    delete canvases[i];

}


file->Close();

delete file;


std::cout << "Valores de mu_values calculados:" << std::endl;

for (float mu : mu_values) {

    std::cout << mu << std::endl;

}


return mu_values;

}

void num_phe_candidatos(const std::vector<float>& mu_values, const std::vector<std::vector<float>>& all_signal_sums) {
const int num_channels = 32;

TFile *file = new TFile("proyecto2_funciones2.root", "UPDATE");

if (!file || file->IsZombie()) {

    std::cerr << "Error: No se pudo abrir el archivo ROOT." << std::endl;

    return;

}


for (int i = 0; i < num_channels; ++i) {

    if (mu_values[i] == 0) {

        std::cout << "Warning: mu_values[" << i << "] = 0 (cannot normalize). Skipping channel." << std::endl;

        continue;

    }


    TH1F* hist = new TH1F(Form("hist_signal_sum_no_gaus_Channel%d", i),

                         Form("Normalized Charge Distribution (Channel %d)", i),

                         100, 0, 5);


    hist->GetXaxis()->SetTitle("Normalized Photoelectron Count");

    hist->GetYaxis()->SetTitle("Entries");


    for (float val : all_signal_sums[i]) {

        hist->Fill(val / mu_values[i]);

    }


    hist->Write("", TObject::kOverwrite);

    delete hist;

}


file->Close();

delete file;

}

void ProcesarAmplitudesMaximas(const std::vector<float>& mu_values) {
    const int num_channels = 32;
    const float bin_width = 8.0;  // ancho de bin en ns
    const int bins_per_record = 1024;

    TH1F* hist_max_pulse_ADC[num_channels];
    TH1F* hist_integrated_charge[num_channels];
    TH1F* hist_photoelectron_count[num_channels];

    // Maximum Pulse Amplitude
    for (int i = 0; i < num_channels; ++i) {
        hist_max_pulse_ADC[i] = new TH1F(Form("MaxPulseAmplitude_ADC_Channel%d", i),
                                         "Maximum Pulse Amplitude per Event",
                                         100, 80, 800);
        hist_max_pulse_ADC[i]->GetXaxis()->SetTitle("ADC");
        hist_max_pulse_ADC[i]->GetYaxis()->SetTitle("Entries");
    }

    // Integrated Pulse Charge
    for (int i = 0; i < num_channels; ++i) {
        hist_integrated_charge[i] = new TH1F(Form("IntegratedPulseCharge_Channel%d", i),
                                            "Integrated Pulse Charge per Event",
                                            100, 0, 5000);
        hist_integrated_charge[i]->GetXaxis()->SetTitle("Charge (ADC*ns)");
        hist_integrated_charge[i]->GetYaxis()->SetTitle("Entries");
    }

    // Estimated Photoelectron Count
    for (int i = 0; i < num_channels; ++i) {
        hist_photoelectron_count[i] = new TH1F(Form("EstimatedPhotoelectronCount_Channel%d", i),
                                              "Estimated Photoelectron Count per Event",
                                              100, 0, 80        );
        hist_photoelectron_count[i]->GetXaxis()->SetTitle("Photoelectron Count");
        hist_photoelectron_count[i]->GetYaxis()->SetTitle("Entries");
    }

    TFile* file = new TFile("proyecto2_funciones2.root", "UPDATE");
    for (int channel = 0; channel < num_channels; ++channel) {
        TTree* tree = (TTree*)file->Get(Form("Channel%d", channel));
        float Data[bins_per_record];
        tree->SetBranchAddress("Data", Data);

        float t_min = 4800.0;
        float t_max = 5000.0;
        int bin_start = static_cast<int>(t_min / bin_width);
        int bin_end = static_cast<int>(t_max / bin_width);

        // Evitar bins fuera de rango
        int bin_start_corrected = std::max(0, bin_start);
        int bin_end_corrected   = std::min(bins_per_record - 1, bin_end);

        for (int entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);

            float baseline = std::accumulate(Data, Data + 300, 0.0f) / 300;
            for (int i = 0; i < bins_per_record; ++i) {
                Data[i] = -(Data[i] - baseline);
            }

            float sum = 0.0;
            float max_amp = -FLT_MAX;
            for (int i = bin_start_corrected; i <= bin_end_corrected; ++i) {
                sum += Data[i];
                if (Data[i] > max_amp) max_amp = Data[i];
            }

            sum *= bin_width;  // multiplicar por ancho de bin para obtener ADC*ns

            if (max_amp > -FLT_MAX) {
                hist_max_pulse_ADC[channel]->Fill(max_amp);
            }
            hist_integrated_charge[channel]->Fill(sum);
            if (mu_values[channel] != 0) {
                hist_photoelectron_count[channel]->Fill(sum / mu_values[channel]);
            } else {
                // Evitar división por cero
                std::cerr << "Warning: mu_values[" << channel << "] = 0, no se puede normalizar." << std::endl;
            }
        }
    }

    for (int i = 0; i < num_channels; ++i) {
        hist_max_pulse_ADC[i]->Write();
        delete hist_max_pulse_ADC[i];
    }
    for (int i = 0; i < num_channels; ++i) {
        hist_integrated_charge[i]->Write();
        delete hist_integrated_charge[i];
    }
    for (int i = 0; i < num_channels; ++i) {
        hist_photoelectron_count[i]->Write();
        delete hist_photoelectron_count[i];
    }

    file->Close();
    delete file;
}


void ProcesarDatosFotoelectrones(const std::vector<float>& mu_values, std::vector<std::vector<float>>& all_signal_sums) {
const int num_channels = 32;
const int bins_per_record = 1024;
const int bin_width = 8.0;

TFile *file = new TFile("proyecto2_funciones2.root", "UPDATE");

if (!file || file->IsZombie()) {

    std::cerr << "Error: No se pudo abrir el archivo ROOT." << std::endl;

    return;

}


// Crear histogramas 2D

TH2F* hist_2D_signal[num_channels];

for (int i = 0; i < num_channels; ++i) {

    hist_2D_signal[i] = new TH2F(Form("hist_2D_signal_Channel%d", i),

                                  Form("ADC vs ADC*ns/mu Channel %d", i),

                                  100, 0, 5,

                                  100, 0, 800);

    hist_2D_signal[i]->GetXaxis()->SetTitle("No. photoelectrons");

    hist_2D_signal[i]->GetYaxis()->SetTitle("Amplitude, ADC");
    
    hist_2D_signal[i]->GetXaxis()->SetNdivisions(5, kTRUE); // el segundo parámetro indica que solo muestra enteros


}


for (int channel = 0; channel < num_channels; ++channel) {

    TTree* tree = (TTree*)file->Get(Form("Channel%d", channel));

    float Data[bins_per_record];

    tree->SetBranchAddress("Data", Data);


    for (int entry = 0; entry < tree->GetEntries(); ++entry) {

        tree->GetEntry(entry);


        float baseline = 0.0;

        int n_baseline = 300;

        for (int i = 0; i < n_baseline; ++i) baseline += Data[i];

        baseline /= n_baseline;


        TH1F *hist_temp = new TH1F("hist_temp", "Histograma Temporal para Baseline", n_baseline, 0, n_baseline);

        for (int i = 0; i < n_baseline; ++i) {

            hist_temp->Fill(Data[i] - baseline);

        }


        float sigma = hist_temp->GetRMS();

        delete hist_temp;


        for (int i = 0; i < bins_per_record; ++i) {

            Data[i] = -(Data[i] - baseline);

        }


        for (int i = 1; i < bins_per_record - 5; ++i) {

            float signal_amplitude = Data[i];

            if (signal_amplitude >= 5 * sigma && signal_amplitude <= 20 * sigma) {

                float sum = 0.0;

                float max_amplitude = 0.0;


                // Calcular suma y encontrar amplitud máxima en la ventana

                for (int j = i - 2; j <= i + 4; ++j) {

                    sum += Data[j];

                    if (Data[j] > max_amplitude) {

                        max_amplitude = Data[j];

                    }

                }
                sum *= bin_width;


                // Calcular el número estimado de fotoelectrones

                float estimated_photoelectrons = sum / mu_values[channel];


                // Llenar el histograma 2D

                hist_2D_signal[channel]->Fill(estimated_photoelectrons, max_amplitude);


                // Almacenar datos

                all_signal_sums[channel].push_back(sum);

                i += 4;

            }

        }

    }

}


// Guardar histogramas 2D

for (int i = 0; i < num_channels; ++i) {

    hist_2D_signal[i]->Write();

    delete hist_2D_signal[i];

}


file->Close();

delete file;

}

int main() {
std::vector<std::vector<float>> all_signal_sums(32);

// Ejecución en orden

CrearAchivoROOT();

GraficarDatosReajustados(all_signal_sums);

GenerarHistograma(all_signal_sums);

std::vector<float> mu_values = gaussianas();

num_phe_candidatos(mu_values, all_signal_sums);

ProcesarAmplitudesMaximas(mu_values);

ProcesarDatosFotoelectrones(mu_values, all_signal_sums);



return 0;

}
