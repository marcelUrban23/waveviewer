#ifndef WAVEDATA_HPP
#define WAVEDATA_HPP

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>
#include <numeric> // std::accumulate


//#define DEBUG_FUNCTION_CALL

class WaveData
{
const std::vector<std::string> m_datatypes {"std::complex<float>",
                                            "std::complex<int32_t>",
                                            "std::complex<int16_t>",
                                            "std::complex<int8_t>",
                                            "float",
                                            "int32_t",
                                            "int16_t",
                                            "int8_t"};
const std::vector<int> m_datatypes_sizeof {8,8,4,2,4,4,2,1};



std::string m_filename;
std::istream *m_input;

std::fstream *m_file;

std::vector<std::complex<float>> m_output;

uint64_t m_filesize_in_bytes = 0;

bool m_swapped = false,
     m_detect_datatype_from_filename = true,
     m_is_file = false;

std::string m_datatype = m_datatypes[0];
uint32_t m_datatype_index = 0;
int32_t  m_byte_offset = 0;

int64_t m_step_width = 0;

public:

double m_samp_rate;


/// @brief Konstruktor
/// @param filename Wenn leer, dann Std::cin
WaveData(std::string filename = {})
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveData::WaveData()" << std::endl;
#endif
    m_file = new std::fstream(filename, std::ios::in | std::ios::binary);
    if ( ! m_file->is_open())
    {
#ifdef DEBUG
        throw std::runtime_error("FEHLER WaveData konnte Datei nicht "
                                 "oeffnen: " + filename);
#endif
        m_input = &std::cin;
    }
    else
    {
        m_input = m_file;

        m_datatype_index = checkFilenameForDatatype(filename);

        m_filesize_in_bytes = m_input->tellg();
        m_input->seekg(0, std::ios::end);
        m_filesize_in_bytes = static_cast<uint64_t>(m_input->tellg())
                            - m_filesize_in_bytes;
        m_input->seekg(0);

        m_filename = filename;
        m_is_file = true;

        m_samp_rate = checkFilenameForSamplerate(filename);
    }

    return;
}

~WaveData()
{
    m_output.resize(0);
    m_file->close();
}


std::vector<std::string> getDatatypes(void){return m_datatypes;}

std::string getFilename(void){return m_filename;}


/// @brief Gibt bestimmten Signalabschnitt zurueck
/// @param begin Startposition relativ zum Anfang
/// @param amount Anzahl inkl. Startposition
/// @return Komplexe Abtastwerte
std::vector<std::complex<float>>
getSamples(uint64_t begin, uint64_t amount)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveData::getsamples()" << std::endl;
#endif
    if( ! m_is_file)
    {
        throw std::runtime_error("FEHLER WaveData::getSamples(): "
                                 "Keine Datei angegeben");
    }
    else
    {
        if(begin <   m_filesize_in_bytes
                   / m_datatypes_sizeof[m_datatype_index])
        {
            m_input->seekg(m_byte_offset
                           + begin * m_datatypes_sizeof[m_datatype_index],
                           std::ios::beg);
        }
    }

    return getNextBlock(amount);
}



/// @brief Gibt die Daten Blockweise aus der Datei zurueck. Der naechste Block
///        wird nach Schrittweite (setNextBlockStepWidth()) gelesen.
/// @param block_leng Anzahl der Samples.
/// @return block_leng complex<float> bzw. weniger wenn Dateiende erreicht.
std::vector<std::complex<float>>
getNextBlock(uint64_t block_leng)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveData::getNextBlock()" << std::endl;
#endif
    if(block_leng != m_output.size())
    {
        m_output.resize(block_leng);
    }

    if((  m_filesize_in_bytes - m_input->tellg())
        / m_datatypes_sizeof[m_datatype_index] < block_leng)
    {
        m_output.resize((  m_filesize_in_bytes - m_input->tellg())
                        / m_datatypes_sizeof[m_datatype_index]);
        block_leng = m_output.size();
    }


    switch(m_datatype_index)
    {
        case 0: // std::complex<float>
        {
             m_input->read(reinterpret_cast<char*>(m_output.data()),
                           block_leng * sizeof(std::complex<float>));

            break;
        }
        case 1: // std::complex<int32_t>
        {
            std::vector<std::complex<int32_t>> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          block_leng * sizeof(std::complex<int32_t>));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](std::complex<int32_t> value)
                           {return std::complex<float>(
                                        static_cast<float>(value.real()),
                                        static_cast<float>(value.imag()));});

            break;
        }
        case 2: // std::complex<int16_t>
        {
            std::vector<std::complex<int16_t>> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          block_leng * sizeof(std::complex<int16_t>));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](std::complex<int16_t> value)
                           {return std::complex<float>(
                                        static_cast<float>(value.real()),
                                        static_cast<float>(value.imag()));});
            break;
        }
        case 3: // std::complex<int8_t>
        {
            std::vector<std::complex<int8_t>> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          tmp.size() * sizeof(std::complex<int8_t>));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](std::complex<int8_t> value)
                           {return std::complex<float>(
                                        static_cast<float>(value.real()),
                                        static_cast<float>(value.imag()));});
            break;
        }
        case 4: // float
        {
            std::vector<float> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          block_leng * sizeof(float));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](float value)
                           {return std::complex<float>(value, 0.0);});
            break;
        }
        case 5: // int32_t
        {
            std::vector<int32_t> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          block_leng * sizeof(int32_t));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](int32_t value)
                           {return std::complex<float>(
                                        static_cast<float>(value), 0.0);});

            break;
        }
        case 6: // int16_t
        {
            std::vector<int16_t> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          block_leng * sizeof(int16_t));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](int16_t value)
                           {return std::complex<float>(
                                        static_cast<float>(value), 0.0);});
            break;
        }
        case 7: // int8
        {
            std::vector<int8_t> tmp(block_leng);
            m_input->read(reinterpret_cast<char*>(tmp.data()),
                          block_leng * sizeof(int8_t));

            std::transform(tmp.begin(), tmp.end(),
                           m_output.begin(),
                           [](int8_t value)
                           {return std::complex<float>(
                                        static_cast<float>(value), 0.0);});
            break;
        }
        default:
        {
            throw std::runtime_error("FEHLER WaveData::getNextBlock(): "
                                     "unbekannter Dateityp");
            break;
        }
    }


    if(m_swapped)
    {
        char* tmp = reinterpret_cast<char*>(m_output.data());

        for(uint64_t w = 0; w <   m_output.size()
                                * sizeof(std::complex<float>); w += 4)
        {
            std::swap(tmp[w    ], tmp[w + 3]);
            std::swap(tmp[w + 1], tmp[w + 2]);
        }
    }

    if(m_step_width)
    {
        m_input->seekg( (m_step_width - block_leng)
                       * m_datatypes_sizeof[m_datatype_index], std::ios::cur);
    }


    return m_output;
}


void setSwapped (bool state){m_swapped = state;}
void setDataType(uint32_t data_type_index){m_datatype_index = data_type_index;}
void setByteOffset(int32_t offset){m_byte_offset = offset;}
int32_t getByteOffset(void){return m_byte_offset;}

/// @brief Setzt die Schrittweite für getNextBlock()
void setNextBlockStepWidth(int64_t step_width){m_step_width = step_width;}


/// @return Anzahl der Abtastwerte in der Datei
std::string getDataType(void)   {return m_datatypes[m_datatype_index];}
uint32_t getDataTypeIndex(void) {return m_datatype_index;}

uint64_t getFileSizeInSamples(void)
{
    if( ! m_is_file)
    {
        throw std::runtime_error("FEHLER WaveData::getFileSizeInSamples: "
                                  "keine Datei geoeffnet");
    }

    return m_filesize_in_bytes / m_datatypes_sizeof[m_datatype_index];
}

/// @brief Setzt die aktuelle Position in der Datei auf die Abtastwert-Position
void setSamplePosition(uint64_t position)
{
    m_input->seekg(position * m_datatypes_sizeof[m_datatype_index],
                   std::ios::beg);
}

/// @return Aktuelle Abtastwert-Position in der Datei
uint64_t
getSamplePosition(void)
{
    if( ! m_is_file)
    {
        throw std::runtime_error("FEHLER WaveData::getFileSizeInSamples: "
                                  "keine Datei geoeffnet");
    }

    return m_input->tellg() * m_datatypes_sizeof[m_datatype_index];
}



/// @brief Gibt Auskunft, von wo die Eingangsdaten gelesen werden.
/// @return true: Eingang ist Datei, false: std::in
bool
getIsInputFile(void)
{
    return m_is_file;
}


/// @brief Prueft das Dateinamenende auf einen Datentyp.
/// @param filename Dateiname
/// @return Index des passenden Datentypeintrags
int
checkFilenameForDatatype(std::string filename)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveData::checkFilenameForDatatype()" << std::endl;
#endif
    if     (filename.find(".32fc") != std::string::npos) return 0;
    else if(filename.find(".32tc") != std::string::npos) return 1;
    else if(filename.find(".16tc") != std::string::npos) return 2;
    else if(filename.find(".8tc")  != std::string::npos) return 3;
    else if(filename.find(".32f")  != std::string::npos) return 4;
    else if(filename.find(".32t")  != std::string::npos) return 5;
    else if(filename.find(".16t")  != std::string::npos) return 6;
    else if(filename.find(".8t")   != std::string::npos) return 7;
    else                                                 return 0;
}

/// @brief Prueft Dateinamen auf enthaltene Abtastrate.
/// @param filename Dateiname
double
checkFilenameForSamplerate(std::string filename)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveData::checkFilenameForSamplerate()" << std::endl;
#endif

    std::string filename_tmp = filename.substr(filename.rfind("/"));
    double samp_rate;

    // KiloHz
    if(filename_tmp.find("k") != std::string::npos)
    {
        std::string tmp;
        int origin = filename_tmp.find("k");
        int pos = origin;

        while(std::isdigit(filename_tmp[--pos]))
        {
            tmp += filename_tmp[pos];
        }
        pos = origin;
        while(std::isdigit(filename_tmp[++pos]))
        {
            tmp += filename_tmp[pos];
        }

        for(int w = pos - origin; w < 4; ++w)
        {
            tmp += "0";
        }
        samp_rate = std::stod(tmp);
    }
    // MegaHz
    else if(filename_tmp.find("m") != std::string::npos)
    {
        std::string tmp;
        int origin = filename_tmp.find("m");
        int pos = origin;

        while(std::isdigit(filename_tmp[--pos]))
        {
            tmp += filename_tmp[pos];
        }
        pos = origin;
        while(std::isdigit(filename_tmp[++pos]))
        {
            tmp += filename_tmp[pos];
        }

        for(int w = pos - origin; w < 7; ++w)
        {
            tmp += "0";
        }
        samp_rate = std::stod(tmp);
    }
    // GigaHz
    else if(filename_tmp.find("g") != std::string::npos)
    {
        std::cerr << filename_tmp << std::endl;
        std::string tmp;
        int origin = filename_tmp.find("g");
        int pos = origin;

        while(std::isdigit(filename_tmp[--pos]))
        {
            tmp += filename_tmp[pos];
        }
        pos = origin;
        while(std::isdigit(filename_tmp[++pos]))
        {
            tmp += filename_tmp[pos];
        }

        for(int w = pos - origin; w < 10; ++w)
        {
            tmp += "0";
        }
        samp_rate = std::stod(tmp);
    }
    // Versucht Abtastrate ohne Buchstaben zu interpretieren
    else
    {
        if(filename_tmp.find(".") != std::string::npos)
        {
            // Zusammenhaengende Ziffern merken und umderehen
            int pos = filename_tmp.find(".");
            std::string tmp, tmp2;
            while(std::isdigit(filename_tmp[--pos]))
            {
               tmp += filename_tmp[pos];
            }
            for(uint64_t w = 1; w <= tmp.size(); ++w)
            {
                tmp2 += tmp[tmp.size() - w];
            }
            samp_rate = std::stod(tmp2);
        }
        else
        {
            samp_rate = 1000000.0;
        }

    }

    return samp_rate;
}


int
estimateMostLikelyDatatype(void)
{
    uint64_t amount_of_samples = 100000;

    if(m_filesize_in_bytes / sizeof(std::complex<float>) < 1000000)
    {
        amount_of_samples = m_filesize_in_bytes / sizeof(std::complex<float>);
    }

    for(uint64_t w = 0; w < m_datatypes.size(); ++w)
    {
        setDataType(w);
        getSamples(0, amount_of_samples);
        std::vector<std::complex<float>> tmp(getSamples(0, amount_of_samples));

        std::complex<float> mean =
         std::accumulate(tmp.begin(), tmp.end(), std::complex<float>(0.0, 0.0));
    }

    return 0;
}

};

#endif // WAVEDATA_HPP
