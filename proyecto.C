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
#include <algorithm>
#include <limits>

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
    const int num_channels = 8;

    // Abrir archivo CSV
    std::ifstream archivo("500LinesExample.csv");
    if (!archivo.is_open()) {
        std::cerr << "No se pudo abrir el archivo." << std::endl;
        return;
    }

    // Crear archivo ROOT
    TFile *file = new TFile("proyecto_funciones.root", "RECREATE");

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
            Bin = std::stoi(columnas[2]);
            Bool = std::stoi(columnas[3]);

            for (int i = 0; i < bins_per_record; ++i) {
                if (!safe_stof(columnas[5 + i], Data[i])) continue;
            }

            // Llenar el TTree correspondiente al canal
            trees[Channel]->Fill();
        }
    }

    // Escribir los árboles en el archivo ROOT
    for (int i = 0; i < num_channels; ++i) {
        trees[i]->Write();
        delete trees[i]; // Liberar memoria
    }

    file->Close();
    delete file; // Liberar memoria
}

void GenerarHistograma(std::vector<std::vector<float>>& all_signal_sums) {
    const int num_channels = 8;
    const int bins_per_record = 1024;
    
    // Abrir el archivo ROOT existente
    TFile *file = new TFile("proyecto_funciones.root", "UPDATE"); // Cambiar a "UPDATE" para agregar datos
    TH1F* hist_signal_sum[num_channels]; 

    for (int i = 0; i < num_channels; ++i) {
        hist_signal_sum[i] = new TH1F(Form("hist_signal_sum_Channel%d", i),
                                       Form("Charge Distribution Channel %d", i),
                                       100, 0, 800);
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
            int n_baseline = 250;
            for (int i = 0; i < n_baseline; ++i) baseline += Data[i];
            baseline /= n_baseline;
            
            // Crear histograma temporal para baseline
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
                    signal_sum.push_back(sum);
                    i += 4;
                }
            }

            // Almacenar la suma de señales para este canal
            all_signal_sums[channel].insert(all_signal_sums[channel].end(), signal_sum.begin(), signal_sum.end());

            for (float sum : signal_sum) {
                hist_signal_sum[channel]->Fill(sum);
            }
        }
    }
    
    // Escribir los histogramas en el archivo ROOT
    for (int i = 0; i < num_channels; ++i) {
        hist_signal_sum[i]->Write(); // Guardar el histograma en el archivo
        delete hist_signal_sum[i]; // Liberar memoria
    }

    file->Close(); // Cerrar el archivo
    delete file; // Liberar memoria
}

std::vector<float> gaussianas() {
    const int num_channels = 8;
    std::vector<float> mu_values(num_channels, 0.0);
    TSpectrum spectrum;

    TFile *file = new TFile("proyecto_funciones.root", "UPDATE"); // Cambiado a UPDATE para guardar los gráficos
    TH1F* hist_signal_sum[num_channels];
    TCanvas* canvases[num_channels];

    for (int i = 0; i < num_channels; ++i) {
        hist_signal_sum[i] = (TH1F*)file->Get(Form("hist_signal_sum_Channel%d", i));
        canvases[i] = new TCanvas(Form("canvas_ch%d", i), Form("Channel %d Peak Analysis", i), 800, 600);
    }

    std::vector<float> target_positions = {100.0};

    for (int i = 0; i < num_channels; ++i) {
        canvases[i]->cd();
        hist_signal_sum[i]->Draw();
        
        // Buscar picos
        int n_peaks = spectrum.Search(hist_signal_sum[i], 1, "nobackground", 0.1);
        hist_signal_sum[i]->Draw(); // Redibujar para mostrar marcadores de picos

        if (n_peaks > 0) {
            Double_t *xpeaks = spectrum.GetPositionX();
            float spe_peak_position = 0.0;

            // Buscar el pico más cercano a 100
            for (float target : target_positions) {
                for (int j = 0; j < n_peaks; ++j) {
                    if (std::abs(xpeaks[j] - target) < 10.0) {
                        spe_peak_position = xpeaks[j];
                        break;
                    }
                }
            }

            if (spe_peak_position != 0.0) {
                float fit_range = 50.0;
                float min_fit = std::max<double>(spe_peak_position - fit_range / 2.0, hist_signal_sum[i]->GetXaxis()->GetXmin());
                float max_fit = std::min<double>(spe_peak_position + fit_range / 2.0, hist_signal_sum[i]->GetXaxis()->GetXmax());
                float height = hist_signal_sum[i]->GetBinContent(hist_signal_sum[i]->FindBin(spe_peak_position));
                float sigma_est = hist_signal_sum[i]->GetRMS();

                TF1 *gaus_fit = new TF1(Form("gaus_spe_%d", i), "gaus", min_fit, max_fit);
                gaus_fit->SetParameters(height, spe_peak_position, sigma_est);
                gaus_fit->SetLineColor(kRed);
                hist_signal_sum[i]->Fit(gaus_fit, "RQ+");
                mu_values[i] = gaus_fit->GetParameter(1);

                // Añadir leyenda
                TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
                leg->AddEntry(hist_signal_sum[i], "Datos", "l");
                leg->AddEntry(gaus_fit, "Ajuste Gaussiano", "l");
                leg->Draw();
                
                // Guardar el canvas en el archivo ROOT
                canvases[i]->Write(Form("peak_analysis_ch%d", i), TObject::kOverwrite);
                
                delete gaus_fit;
                delete leg;
            } else {
                std::cerr << "Advertencia: No se encontró un pico adecuado en el canal " << i << std::endl;
                hist_signal_sum[i]->SetTitle(Form("Channel %d - No peak found near 100", i));
                canvases[i]->Write(Form("peak_analysis_ch%d", i), TObject::kOverwrite);
            }
        } else {
            std::cerr << "Advertencia: No se encontraron picos en el canal " << i << std::endl;
            hist_signal_sum[i]->SetTitle(Form("Channel %d - No peaks found", i));
            canvases[i]->Write(Form("peak_analysis_ch%d", i), TObject::kOverwrite);
        }
    }

    // Limpiar
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
    const int num_channels = 8;

    TFile *file = new TFile("proyecto_funciones.root", "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: No se pudo abrir el archivo ROOT." << std::endl;
        return;
    }

    for (int i = 0; i < num_channels; ++i) {
        if (mu_values[i] == 0) {
            std::cout << "Warning: mu_values[" << i << "] = 0 (cannot normalize). Skipping channel." << std::endl;
            continue;
        }

        std::cout << "Channel " << i << " - signal sums count: " << all_signal_sums[i].size() << std::endl;

        // Calculate normalized range
        float min_val = std::numeric_limits<float>::max();
        float max_val = std::numeric_limits<float>::lowest();
        
        for (float val : all_signal_sums[i]) {
            float normalized_charge = val / mu_values[i];
            if (normalized_charge < min_val) min_val = normalized_charge;
            if (normalized_charge > max_val) max_val = normalized_charge;
        }
        std::cout << "Channel " << i << ": normalized range = [" << min_val << ", " << max_val << "]" << std::endl;

        // Create normalized histogram
        TH1F* hist = new TH1F(Form("hist_signal_sum_no_gaus_Channel%d", i),
                             Form("Normalized Charge Distribution (Channel %d)", i),
                             100, 0, 5);

        hist->GetXaxis()->SetTitle("Normalized Photoelectron Count (sum/mu)");
        hist->GetYaxis()->SetTitle("Entries");

        // Fill histogram with ALL normalized values (no filtering)
        for (float val : all_signal_sums[i]) {
            hist->Fill(val / mu_values[i]);
        }

        hist->Write("", TObject::kOverwrite);
        delete hist;
    }

    file->Close();
    delete file;
}

void CalcularYGuardarAmplitudMaxima(float Data[], int Channel, TH1F* hist_amplitude_max_ADC[], TH1F* hist_amplitude_max[], TH1F* hist_amplitude_max_sum[], const std::vector<float>& mu_values) {
    const float bin_width = 8.0; // Asegúrate de que esto coincida con tu definición
    float t_min = 4800.0;
    float t_max = 5100.0;
    int bin_start = static_cast<int>(t_min / bin_width);
    int bin_end = static_cast<int>(t_max / bin_width);

    if (bin_start < 0) bin_start = 0;
    if (bin_end > 1023) bin_end = 1023; // Asegúrate de que este valor coincida con bins_per_record - 1

    // Inicializar la suma y el valor máximo
    float sum = 0.0;
    float max_amplitude = -std::numeric_limits<float>::infinity(); // Inicializa con el menor valor posible

    for (int i = bin_start; i <= bin_end; ++i) {
        sum += Data[i];
        if (Data[i] > max_amplitude) {
            max_amplitude = Data[i]; // Actualiza el máximo si se encuentra un nuevo valor mayor
        }
    }

    // Llenar el histograma hist_amplitude_max_ADC solo con el valor máximo
    if (max_amplitude > -std::numeric_limits<float>::infinity()) { // Verifica que se haya encontrado un máximo válido
        hist_amplitude_max_ADC[Channel]->Fill(max_amplitude);
    }

    hist_amplitude_max[Channel]->Fill(sum / mu_values[Channel]);
    hist_amplitude_max_sum[Channel]->Fill(sum);
}

void ProcesarAmplitudesMaximas(const std::vector<float>& mu_values) {
    const int num_channels = 8;
    const float bin_width = 8.0;
    const int bins_per_record = 1024;

    // Crear histogramas
    TH1F* hist_amplitude_max_ADC[num_channels];
    TH1F* hist_amplitude_max[num_channels];
    TH1F* hist_amplitude_max_sum[num_channels];

    // Sistema de nombres para agrupar en TBrowser:
    // Usar prefijos consistentes para agrupar por tipo
    
    // 1. Primer grupo: Todos los ADC (A_*)
    for (int i = 0; i < num_channels; ++i) {
        hist_amplitude_max_ADC[i] = new TH1F(Form("A_ADC_Channel%d", i), 
                                           "Max Amplitude ADC", 100, 0, 200);
        hist_amplitude_max_ADC[i]->GetXaxis()->SetTitle("ADC");
        hist_amplitude_max_ADC[i]->GetYaxis()->SetTitle("Entries");
    }

    // 2. Segundo grupo: Todos los Max (M_*)
    for (int i = 0; i < num_channels; ++i) {
        hist_amplitude_max[i] = new TH1F(Form("M_Max_Channel%d", i),
                                      "Charge Distribution Max", 100, 2, 10);
        hist_amplitude_max[i]->GetXaxis()->SetTitle("Photoelectron Count");
        hist_amplitude_max[i]->GetYaxis()->SetTitle("Entries");
    }

    // 3. Tercer grupo: Todos los Sum (S_*)
    for (int i = 0; i < num_channels; ++i) {
        hist_amplitude_max_sum[i] = new TH1F(Form("S_Sum_Channel%d", i),
                                          "Sum of Amplitudes", 100, 0, 800);
        hist_amplitude_max_sum[i]->GetXaxis()->SetTitle("Charge ADC*ns");
        hist_amplitude_max_sum[i]->GetYaxis()->SetTitle("Entries");
    }

    TFile* file = new TFile("proyecto_funciones.root", "UPDATE");
    for (int channel = 0; channel < num_channels; ++channel) {
        TTree* tree = (TTree*)file->Get(Form("Channel%d", channel));
        float Data[bins_per_record];
        tree->SetBranchAddress("Data", Data);

        float t_min = 4800.0;
        float t_max = 5100.0;
        int bin_start = static_cast<int>(t_min / bin_width);
        int bin_end = static_cast<int>(t_max / bin_width);

        for (int entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);
            
            float baseline = std::accumulate(Data, Data + 250, 0.0f) / 250;
            for (int i = 0; i < bins_per_record; ++i) {
                Data[i] = -(Data[i] - baseline);
            }

            float sum = 0.0;
            float max_amp = -FLT_MAX;
            for (int i = bin_start; i <= bin_end; ++i) {
                sum += Data[i];
                if (Data[i] > max_amp) max_amp = Data[i];
            }

            if (max_amp > -FLT_MAX) {
                hist_amplitude_max_ADC[channel]->Fill(max_amp);
            }
            hist_amplitude_max[channel]->Fill(sum / mu_values[channel]);
            hist_amplitude_max_sum[channel]->Fill(sum);
        }
    }

    // Guardar en el orden deseado
    for (int i = 0; i < num_channels; ++i) {
        hist_amplitude_max_ADC[i]->Write();
    }
    for (int i = 0; i < num_channels; ++i) {
        hist_amplitude_max[i]->Write();
    }
    for (int i = 0; i < num_channels; ++i) {
        hist_amplitude_max_sum[i]->Write();
    }

    // Limpiar memoria
    for (int i = 0; i < num_channels; ++i) {
        delete hist_amplitude_max_ADC[i];
        delete hist_amplitude_max[i];
        delete hist_amplitude_max_sum[i];
    }
    file->Close();
    delete file;
}



int main() {
    std::vector<std::vector<float>> all_signal_sums(8); // 8 canales
    // Ejecución en orden
    CrearAchivoROOT();
    GenerarHistograma(all_signal_sums);
    std::vector<float> mu_values = gaussianas();
    num_phe_candidatos(mu_values, all_signal_sums);
    ProcesarAmplitudesMaximas(mu_values);
    return 0;
}