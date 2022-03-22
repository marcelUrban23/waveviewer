#ifndef libdsss_HPP
#define libdsss_HPP

#include <vector>
#include <numeric>
#include <complex>
#include <iostream>
#include <algorithm>

#include "fft2.hpp"
#include "fftwindows.hpp"


//#define DEBUG_FUNCTION_CALL

class Dsss
{
private:

/* Vektor mit spreitz Sequenzen. */
std::vector<std::vector<int>>                 m_spread_sequences,
                                              m_lsr_sequences;

/* Vektor mit fft transformierte spreitz-Sequenzen. */
std::vector<std::vector<std::complex<float>>> m_fft_sequences,
                                              m_xcorr_result;


double   m_peak_to_average_ratio;


uint64_t m_fft_exp,
         m_fft_leng;

FFT2 *m_fft2;

FFTWindow *m_fft_window;

public:

Dsss(uint32_t fft_exp)
{
    m_spread_sequences.reserve(100);
    m_lsr_sequences.reserve(100);
    m_fft_sequences.reserve(100);
    m_xcorr_result.reserve(100);

    /* Default Erkennungs-Schwelle setzen */
    m_peak_to_average_ratio = 7.0;

    m_fft_exp = fft_exp;
    m_fft_leng = 1 << m_fft_exp;


    m_fft2 = new FFT2(fft_exp);

    m_fft_window = new FFTWindow();

    return;
}
~Dsss()
{
    m_spread_sequences.resize(0);
    m_lsr_sequences.resize(0);
    m_fft_sequences.resize(0);
    m_xcorr_result.resize(0);

    delete m_fft2;
}




/// @brief Fuegt einen zu pruefenden Chip hinzu.
/// @param Polynom des LSR
/// @param Vorbelegung des LSR
/// @param lsr_sequence_leng Laenge der
/// @param bit_repeat Anzahl der Wiederholungen jeder Polaritaet
void
addSequence(std::vector<int>  lsr_polynom,
            std::vector<int>  lsr_seed,
            uint32_t          lsr_sequence_leng,
            uint32_t          bit_repeat)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::addSequence()\n";
#endif

    /* Schieberegisterfolge mit Vorbelegung erzeugen. */
    m_lsr_sequences.push_back(genLSRSequence(lsr_polynom,
                                             lsr_seed,
                                             lsr_sequence_leng));

    m_spread_sequences.push_back(genSpreadSequence(m_lsr_sequences.back(),
                                                   bit_repeat));

    std::vector<std::complex<float>> tmp(m_fft_leng);

    for(uint64_t w = 0; w < m_spread_sequences.back().size(); ++w)
    {
        tmp[w].real(m_spread_sequences.back()[w]);
    }


    /* Foutier Transformierte der generierten Sequenz bilden. */
    m_fft_sequences.push_back(m_fft2->fft(tmp));


    /* neuen Ausgangsslot fuer die Berechnungsergebnisse hinzufuegen. */
    std::vector<std::complex<float>> new_xcorr_result(m_fft_leng);
    m_xcorr_result.push_back(new_xcorr_result);


    return;
}



/// @brief Korreliert alle hinzugefuegten Sequenzen mit dem Eingangsvektor.
std::vector<int>
correlateSequences(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::correlateSequences()\n";
#endif
    std::vector<std::complex<float>> input_fft(m_fft2->fft(input));

    std::vector<int>     result(m_fft_sequences.size());

    for(uint64_t x = 0; x < m_fft_sequences.size(); ++x)
    {
        /* Komplex konjugierte Mukltiplikation mit Sequenz durchfueheren. */
        for(uint64_t w = 0; w < m_fft_leng; ++w)
        {
            m_xcorr_result[x][w] = std::conj(input_fft[w])
                                 * m_fft_sequences[x][w];
        }

        /* Ergbnis vom Frequenzbereich in den Zeitbereich transformieren. */
        m_xcorr_result[x] = m_fft2->fft(m_xcorr_result[x], 1);


        //fft_input, m_fft_leng, m_xcorr_result[x].begin()

    }
    return result;
}



std::vector<int>
genLSRSequence(const std::vector<int> &lsr_polynom,
               std::vector<int>        lsr_seed,
               uint32_t                lsr_sequence_leng)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::genLSRSequence()\n";
#endif

    uint32_t polynom_degree = 0,
             w              = 0,
             x              = 0;


    /* Polynomgrad bestimmen. */
    polynom_degree = *std::max_element(lsr_polynom.begin(), lsr_polynom.end());

    /* Speicher mit 0 initialisieren. */
    std::vector<int> lsr_sequence(lsr_sequence_leng, 0);

    /* Vorbelegung setzen. */
    std::copy(lsr_seed.begin(),
              lsr_seed.begin() + lsr_seed.size(),
              lsr_sequence.begin());


    /* LSR Folge bis zur gesetzten Laenge durchtakten. */
    for(w = 0; w < (lsr_sequence_leng - polynom_degree); ++w)
    {
        /* Alle Abgriffe gemaeß Exponenten in den .real() schreiben. */
        for(x = 0; x < lsr_polynom.size(); ++x)
        {
            lsr_sequence[polynom_degree + w] = lsr_sequence[polynom_degree + w]
                                             ^ lsr_sequence[w + lsr_polynom[x]];
        }
    }


    return lsr_sequence;
}

/* Folgend werden alle Werte von 0 auf -1 geaendert: */
/* z.B. 1,1,0,0,1,0 --> +1,+1,-1,-1,+1,-1 */
/* Expandiert jedes Bit einer LSR Sequenz auf die Anzahl bit_repeat. */
/* z.B. fuer bit_repeat = 2 wird jedes Bit verdoppelt: */
/* +1,+1,-1,-1,+1,-1 --> +1,+1,+1,+1,-1,-1,-1,-1,+1,+1,-1,-1 */
std::vector<int>
genSpreadSequence(const std::vector<int> &sequence, uint64_t bit_repeat)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::genSpreadSequence()\n";
#endif
    std::vector<int> result(sequence.size() * bit_repeat);

    /* Gehe alle Registerwerte durch. */
    for(uint64_t w = 0; w < sequence.size(); ++w)
    {
        /* Vervielfache jeden Registerwert. */
        for(uint64_t x = 0; x < bit_repeat; ++x)
        {
            sequence[w] == 0 ? (result[w * bit_repeat + x] = -1)
                             : (result[w * bit_repeat + x] =  1);
        }
    }


    return result;
}



/******************************************************************************/
static int32_t
findDominantPeak(const std::vector<double> &input,
                 uint32_t            guard_interval,
                 double             *peak_to_side_ratio)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::findDominantPeak()\n";
#endif


    int64_t  space_to_peak = 5,  /* Sicherheitsabstand */
             area_to_peak  = 25;  /* Beachtete Werte nahe der Spitze */



    /* Randeffekte ausblenden indem ein Schutzabstand gehalten wird. */
    std::vector<double>::const_iterator result =
                           std::max_element(input.begin() + guard_interval,
                                            input.end()   - guard_interval);

    double peak = *result;

    int32_t peak_pos = std::distance(input.begin(), result);

    double side_average = 0.0;
    /* Aufaddieren der Betraege links und rechts der Spitze. */
    for(int64_t w = 0; w < area_to_peak; w++)
    {
        side_average += (  input[peak_pos - space_to_peak - w]
                         + input[peak_pos + space_to_peak + w]);
    }

    /* Aritmetisches Mittel bilden.*/
    side_average /= (2 * area_to_peak);

    /* Spitze-zu-Umgebungs Verhaeltnis errechnen. */
    peak_to_side_ratio[0] = peak / side_average;


    return peak_pos;
}



void setPeakToAverageRatio(double ratio){m_peak_to_average_ratio = ratio;}
double getPeakToAverageRatio(void){return m_peak_to_average_ratio;}
std::vector<int> getSequence(int index){return m_lsr_sequences[index];}
std::vector<std::complex<float>> getFFTSequence(int index)
{
    return m_fft_sequences[index];
}


int getNumOfSequences(void){
    return m_spread_sequences.size();
}



std::vector<std::complex<float>>
getXcorrResult(int index)
{
    return m_xcorr_result[index];
}




/*************************/
/* Allgemeine Funktionen */
/*************************/






/// @brief berechnet die Kreuzkorrelierte zweier Signale
/// @param input Abtastwerte Signal 1
/// @param input Abtastwerte Signal 2
/// @return Kreuzkorrelierte der Eingangsdaten
template <typename T = std::complex<float>>
std::vector<T>
getCrossCorr(const std::vector<T> &input,
             const std::vector<T> &sequence)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getCrossCorr(1)\n";
#endif
    std::vector<T> result(m_fft2->fft(input));

    /* Komplex konjugierte Multiplikation */
    std::transform(result.begin(),
                   result.end(),
                   m_fft2->fft(sequence).begin(),
                   result.begin(),
                   [](std::complex<float> a,
                      std::complex<float> b)
                      {return a * std::conj(b);});

    result = m_fft2->fft(result, 1);


    return result;
}


/// @brief
std::vector<std::complex<float>>
getCrossCorr(const std::vector<std::complex<float>> &input,
             const std::vector<int>                 &sequence)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getCrossCorr(2)\n";
#endif
    std::vector<std::complex<float>> tmp(input.size(), std::complex<float>(0.0, 0.0));
    for(uint64_t w = 0; w < sequence.size(); ++w)
    {
        tmp[w] = std::complex<float>(sequence[w], 0.0);
    }

    std::vector<std::complex<float>> result(m_fft2->fft(input));

    /* Komplex konjugierte Multiplikation */
    std::transform(result.begin(),
                   result.end(),
                   m_fft2->fft(tmp).begin(),
                   result.begin(),
                   [](std::complex<float> a,
                      std::complex<float> b)
                      {return a * std::conj(b);});

    result = m_fft2->fft(result, 1);


    return result;
}

template <typename T = std::complex<float>>
static double
getVariance(const std::vector<T> &input)
{
    double sum = 0.0;
    for(auto val : input)
    {
        // Quadrat des Betrages (std::pow(std::abs, 2))
        sum += static_cast<double>(std::real(val * std::conj(val)));
    }

    return sum / input.size();
}

static uint64_t
getNextPow(uint64_t leng)
{
    uint64_t exp = 0;
    while(1 << ++exp < leng);
    return exp;
}


/// @brief
std::vector<double>
getNormalizedCrossCorr(const std::vector<std::complex<float>> &input,
                       const std::vector<std::complex<float>> &sequence)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdss::getNormalizedCrossCorr(1)" << std::endl;
#endif

    std::vector<std::complex<float>> in_fft(m_fft2->fft(input)),
                                     seq_fft(input.size(),
                                             std::complex<float>(0.0, 0.0)),
                                     result(input.size());

    std::copy(sequence.begin(), sequence.end(), seq_fft.begin());
    seq_fft = m_fft2->fft(seq_fft);

    float a = 0.0;
    for(uint64_t w = 0; w < seq_fft.size(); ++w)
    {
        a += std::norm(seq_fft[w]);
    }
    std::transform(result.begin(), result.end(), result.begin(),
                   [a](std::complex<float> value)
                   {return value / a;});


    /* Komplex konjugierte Multiplikation */
    std::transform(in_fft.begin(), in_fft.end(),
                   seq_fft.begin(),
                   result.begin(),
                   [](std::complex<float> a,
                      std::complex<float> b)
                      {return a * std::conj(b);});

    result = m_fft2->fft(result, 1);

    /* Normiert auf die lokale Energy. */
    for(uint64_t w = 0; w < result.size() - sequence.size(); ++w)
    {
        double local_energy = 0.0;
        for(uint64_t x = 0; x < sequence.size(); ++x)
        {
            local_energy += std::norm(result[w + x]);
        }
        result[w] /= local_energy;
    }

    for(uint64_t w = result.size() - sequence.size(); w < result.size(); ++w)
    {
        result[w] = std::complex<float>(0.0, 0.0);
    }

//    double a = 0.0, b = 0.0;
//    for(uint64_t w = 0; w < input.size(); ++w)
//    {
//        a += std::norm(in_fft[w]);
//        b += std::norm(seq_fft[w]);
//    }

//    float div = std::sqrt(a * b) / input.size();

//    std::vector<double> out(input.size());
//    std::transform(result.begin(), result.end(), out.begin(),
//                   [div](std::complex<float> value) -> double
//                   {return std::abs(value / div);});


    std::vector<double> out(input.size());
    std::transform(result.begin(), result.end(), out.begin(),
                   [](std::complex<float> value) -> double
                   {return std::abs(value);});


    return out;
}



std::vector<double>
getNormalizedCrossCorr(const std::vector<std::complex<float>> &input,
                       const std::vector<int>                 &sequence)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getNormalizedCrossCorr(2)" << std::endl;
#endif
    if(input.size() < sequence.size()) return {};

    std::vector<std::complex<float>> seq(sequence.size(),
                                         std::complex<float>(0.0, 0.0));
    for(uint64_t w = 0; w < sequence.size(); ++w)
    {
        seq[w].real(sequence[w]);
    }


    return getNormalizedCrossCorr(input, seq);
}




/// @brief berechnet die Autokorrelierte der Eingangsdaten
///        ->fft->konjKomplMult->fftInvers
/// @param input Abtastwerte
/// @return Autokorrelierte der Eingangsdaten
std::vector<std::complex<float>>
getAutoCorr(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getAutoCorr()" << std::endl;
#endif
    if(input.size() != m_fft_leng)
    {
        throw std::runtime_error("FEHLER libdsss::getAutoCorr: "
                                 "input.size() != m_fft_leng");
    }

    return m_fft2->fft(getConj(m_fft2->fft(input)), 1);
}




/// @brief Berechnet den Absolutwert
/// @param input Abtastwerte
/// @return Absolutwerte der Eingangsdaten
template <typename T = double>
static std::vector<std::vector<T>> inline
getAbs(const std::vector<std::vector<std::complex<float>>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getAbs(1)" << std::endl;
#endif
    std::vector<std::vector<T>> output(input.size());
    for(uint64_t w = 0; w < input.size(); ++w)
    {
        output[w] = getAbs(input[w]);
    }

    return output;
}
/// @brief Berechnet den Absolutwert
/// @param input Abtastwerte
/// @return Absolutwerte der Eingangsdaten
template <typename T = double>
static std::vector<T> inline
getAbs(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getAbs(2)" << std::endl;
#endif
    std::vector<T> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                   [](std::complex<float> c) -> T
                   {return static_cast<T>(std::abs(c));});


    return output;
}



/// @brief Addiert to_add auf base zellenweise.
template <typename T = std::complex<float>>
static inline void
addToVector(std::vector<T> &base, const std::vector<T> &to_add)
{
    if(base.size() != to_add.size()){base = to_add; return;}

    std::transform(base.begin(), base.end(), to_add.begin(),
                   base.begin(), std::plus<T>());

    return;
}
template <typename T = std::complex<float>>
static inline void
addToVector(      std::vector<std::vector<T>> &base,
            const std::vector<std::vector<T>> &to_add)
{
    if(base.size() != to_add.size()){base = to_add; return;}

    for(uint64_t w = 0; w < base.size(); ++w)
    {
        addToVector(base[w], to_add[w]);
    }

    return;
}


/// @brief Berechnet das Leistungsdichtesprektrum der Eingagsdaten, welche
///        vorher potenziert werden koennen
/// @param input Abtastwerte
/// @param pow Potenz der Eingangsdaten
/// @return Leistungsdichtespektrum
std::vector<double> inline
getPSD(const std::vector<std::complex<float>> &input, uint64_t pow = 1,
       int window_index = 0)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdss::getPSD(2)" << std::endl;
#endif
    return getAbs(getConj(m_fft2->fft(
                                  getWindowedSignal(
                                  getNPow<std::complex<float>>(input, pow),
                                  window_index))));
}


std::vector<std::complex<float>>
getWindowedSignal(const std::vector<std::complex<float>> &input, int index)
{
    if(index == 0)
    {
        return m_fft_window->blackmanWindow(input);
    }
    else if(index == 1)
    {
        return m_fft_window->flattopWindow(input);
    }
    else if(index == 2)
    {
        return m_fft_window->hammingWindow(input);
    }
    else
    {
        return m_fft_window->vonHannWindow(input);
    }
}


/// @brief Gibt den Betrag des fft gefensterten Betrages des Signals zurueck.
/// @param input
std::vector<double>
getAMDetect(const std::vector<std::complex<float>> &input, int window_index = 0)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdss::getAMDetect()" << std::endl;
#endif
    return getAbs<double>(
                m_fft2->fft(
                    getAbs<std::complex<float>>(
                        getWindowedSignal(input,
                                                        window_index))));
}



/// @brief Enspreizt die Eingangsdaten anhand des Chips forlaufend
/// @param input zu entspreizende Abtastwerte
/// @param chip Spreizfolge
/// @return entspreizte Abtastwerte
std::vector<std::complex<float>>
despread(const std::vector<std::complex<float>> &input,
         const std::vector<int> &chip)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdss::despread()" << std::endl;
#endif
    double ratio = 0.0;

    uint64_t pos = findDominantPeak(getAbs(getCrossCorr(input, chip)),
                                    static_cast<uint32_t>(chip.size()),
                                    &ratio);

    std::cerr << ratio << std::endl;

    std::vector<std::complex<float>> output(input.size());

    pos = (chip.size() - (pos % chip.size())) % chip.size();

    /* oldschool: Es wird staendig geprueft, ob die letzte Chipposition *
     * erreicht und wieder vorn gezaehlt werden muss. */
    for(uint64_t w = 0; w < input.size(); ++w)
    {
        chip[pos] < 0 ? (output[w] = -input[w]) : (output[w] = input[w]);

        if(++pos > (chip.size() - 1)) pos = 0;
    }

    return output;
}



std::vector<int>
estimateAllChipLeng(const std::vector<double> &input, double threshold = 5.0)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::estimateBestChipLeng()" << std::endl;
#endif

    double average = std::accumulate(input.begin(), input.end(), 0.0)
                   / input.size();

    uint64_t x = 0;

    std::vector<int> tmp_0;
    for(uint64_t w = 10; w < (input.size() / 2); ++w)
    {
        if(   (input[w] > input[w - 1])
           && (input[w] > input[w + 1])
           && (input[w] > (average * threshold)))
        {
            tmp_0.push_back(w - x);

            x = w;
        }
    }

    if(tmp_0.empty())
    {
        std::cerr << "keine chiplaenge gefunden\n";

        return {};
    }

    /* Aufsteigend sortieren. */
    std::sort(tmp_0.begin(),
              tmp_0.end(),
              std::less<int>());

    std::vector<int> tmp_1;
    /* Alle unterschiedlichen Laengen jeweils einmal merken. */
    for(auto w : tmp_0)
    {
        int is_new = 1;
        for(auto x : tmp_1)
        {
            if(w == x)
            {
                is_new = 0;
            }
        }
        if(   is_new
           && (w > 20))
        {
            tmp_1.push_back(w);
        }
    }

    /* Vorkommen der chiplaengen ermitteln. */
    std::vector<int> tmp_2(tmp_1.size());
    for(uint64_t w = 0; w < tmp_0.size(); ++w)
    {
        for(uint64_t x = 0; x < tmp_1.size(); ++x)
        {
            /* Nur kleinste gemeinsame Teiler ermitteln. */
            if(0 == (tmp_0[w] % tmp_1[x]))
            {
                tmp_2[x]++;
                break;
            }
        }
    }

//#ifdef DEBUG
    std::cerr << "estimateBestChipLeng()\n";
    for(uint64_t w = 0; w < tmp_1.size(); ++w)
    {
        std::cerr << "chip_leng: " << tmp_1[w] << ", hits: " <<
                     tmp_2[w] << std::endl;
    }
//#endif

    std::vector<int> result;

    for(uint64_t w = 0; w < tmp_1.size(); ++w)
    {
        if(tmp_2[w])
        {
            result.push_back(tmp_1[w]);
        }
    }


    /* ermittelte chiplaenge zurueckgeben. */
    return result;
}



/// @brief versucht eine Chiplaenge zu bestimmen
/// @param input Autokorrelation
/// @param threshold >= peak / average
/// @return beste Chiplaenge
int
estimateBestChipLeng(const std::vector<double> &input, double threshold = 5)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::estimateBestChipLeng()" << std::endl;
#endif

    double average = std::accumulate(input.begin(), input.end(), 0.0)
                   / input.size();

    uint64_t x = 0;

    std::vector<uint64_t> tmp_0,
                          tmp_1;

    tmp_0.reserve(input.size());
    tmp_1.reserve(input.size());

    for(uint64_t w = 10; w < (input.size() / 2); ++w)
    {
        if(   (input[w] > input[w - 1])
           && (input[w] > input[w + 1])
           && (input[w] > (average * threshold)))
        {
            tmp_0.push_back(w - x);

            x = w;
        }
    }

    if(tmp_0.empty())
    {
        std::cerr << "keine chiplaenge gefunden\n";

        return -1;
    }

    /* Aufsteigend sortieren. */
    std::sort(tmp_0.begin(),
              tmp_0.end(),
              std::less<uint64_t>());

    /* Alle unterschiedlichen Laengen jeweils einmal merken. */
    for(auto w : tmp_0)
    {
        int is_new = 1;
        for(auto x : tmp_1)
        {
            if(w == x)
            {
                is_new = 0;
            }
        }
        if(   is_new
           && (w > 20))
        {
            tmp_1.push_back(w);
        }
    }

    /* Vorkommen der chiplaengen ermitteln. */
    std::vector<uint64_t> tmp_2(tmp_1.size());


    for(uint64_t w = 0; w < tmp_0.size(); ++w)
    {
        for(uint64_t x = 0; x < tmp_1.size(); ++x)
        {
            /* Nur kleinste gemeinsame Teiler ermitteln. */
            if(0 == (tmp_0[w] % tmp_1[x]))
            {
                tmp_2[x]++;
                break;
            }
        }
    }

#ifdef DEBUG
    std::cerr << "estimateBestChipLeng()\n";
    for(uint64_t w = 0; w < tmp_1.size(); ++w)
    {
        std::cerr << "chip_leng: " << tmp_1[w] << ", hits: " <<
                     tmp_2[w] << std::endl;
    }
#endif

    /* Beste Trefferanzehl bestimmen. */
    std::vector<uint64_t>::iterator best = std::max_element(tmp_2.begin(),
                                                            tmp_2.end());
    /* Position besten Trefferanzahl bestimmen. */
    uint64_t dist = std::distance(tmp_2.begin(), best);


    /* ermittelte chiplaenge zurueckgeben. */
    return tmp_1[dist];
}



/// @brief versucht die Polaritaet (1, -1) anhand eine Chiplaenge
///        zu bestimmen
/// @param input Abtastwerte
/// @param pattern_leng angenommene Chiplaenge
/// @return Chip
std::vector<int>
estimateChipPattern(const std::vector<std::complex<float>> &input, int pattern_leng)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::estimateChipPattern()" << std::endl;
#endif


    if(pattern_leng < 20)
    {
        std::cerr << "FEHLER estimateChipPattern(): scheiss chiplaenge\n";

        return {};
    }

    uint64_t chip_cnt = input.size() / pattern_leng;
    std::vector<std::complex<float>> polars(pattern_leng, 0.0);
    for(uint64_t w = 0; w < chip_cnt; ++w)
    {
        for(int x = 0; x < pattern_leng; ++x)
        {
            polars[x] += input[w * pattern_leng + x];
        }
    }

    /* Durchschnittliche Abweichung zu (0.0, 0.0) ermitteln. */
    std::complex<float> mass_point = std::accumulate(polars.begin(),
                                                     polars.end(),
                                                  std::complex<float>(0.0,0.0));
    mass_point /= polars.size();
    /* Orientierung bestimmen. */
    float delta_angel = 0.0;

    /* Durchschnittliche Abweichung zum (0.0, 0.0) abziehen. */
    for(int w = 0; w < pattern_leng; ++w)
    {
        polars[w] -= mass_point;
        delta_angel += std::arg(std::complex<float>(std::abs(polars[w].real()),
                                                   std::abs(polars[w].imag())));

    }

    delta_angel /= polars.size();

    std::vector<int> result(pattern_leng, 0);
    for(int w = 0; w < pattern_leng; ++w)
    {
        (std::polar(std::abs(polars[w]), static_cast<float>(std::arg(polars[w]) - delta_angel))).real()
         < 0.0 ? (result[w] = -1) : (result[w] = 1);
    }


    return result;
}



/// @brief Exponiert jeden Wert pow-fach
/// @param input Eingangsvektor
/// @param pow Exponent
/// @param pow Exponent um den jeder Wert potentiert wird
template <typename T = std::complex<float>>
static std::vector<T> inline
getNPow(const std::vector<T> &input, int pow = 2)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getPow(2)" << std::endl;
#endif
    std::vector<T> result(input.size());
    std::transform(input.begin(), input.end(), result.begin(),
                   [pow](T value) -> T
                   {return std::pow(value, pow);});

    return result;
}



/// @brief Berechnet fuer jeden Wert des Vektors: 10.0 * std::log10(Value).
/// @param input Eingangsvektor
/// @return Ausgangsvektor double
template <typename T = double>
static std::vector<T>
getLog10(const std::vector<float> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getLog10()" << std::endl;
#endif
    std::vector<T> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                   [](float val) -> T
                   {return 10.0 * std::log10(val);});


    return output;
}
/// @brief Berechnet fuer jeden Wert des Vektors: 10.0 * std::log10(Value).
/// @param input Eingangsvektor
/// @return Ausgangsvektor double
template <typename T = double>
static std::vector<T>
getLog10(const std::vector<double> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getLog10()" << std::endl;
#endif
    std::vector<T> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                   [](double val) -> T
                     {return 10.0 * std::log10(val);});;
    return output;
}
/// @brief Berechnet fuer jeden Wert des Vektors: 10.0 * std::log10(Value).
/// @param input Eingangsvektor
/// @return Ausgangsvektor double
template <typename T = double>
static std::vector<T>
getLog10(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getLog10()" << std::endl;
#endif
    std::vector<T> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(),
                   [](std::complex<float> val) -> T
                   {return 10.0 * std::log10(std::abs(val));});


    return output;
}
///// @brief Berechnet fuer jeden Wert des Vektors: 10.0 * std::log10(Value).
///// @param input Eingangsvektor 2D
///// @return Ausgangsvektor 2D double
static std::vector<std::vector<double>>
getLog10(const std::vector<std::vector<std::complex<float>>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getLog10(2)" << std::endl;
#endif
    std::vector<std::vector<double>> output(input.size());
    for(uint64_t w = 0; w < input.size(); ++w)
    {
        output[w] = getLog10(input[w]);
    }

    return output;
}




/// @brief Normiert alle Daten auf die Energie
/// @param input Eingangsvektor
/// @return Normierte Eingangswerte
template <typename T = std::complex<float>>
std::vector<T>
getNormalizedSignal(const std::vector<T> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getNormalizedSignal()" << std::endl;
#endif
    std::vector<T> output(input.size());

    float energy = 0.0;
    for(auto value : input)
    {
        energy += std::norm(value);
    }
    energy = std::sqrt(energy / input.size());

    std::transform(input.begin(), input.end(), output.begin(),
                   [energy](T value)
                   {return value / energy;});

    return output;
}



/// @brief Bestimmt den Vektor mit Traegerfrequenzen des Eingangssignals
/// @param input Eingangssignal
/// @param n Anzahl der Phasenzustaende (wird als Exponent verwendet)
/// @return Leistungsdichtespektrum (PSD)
std::vector<double> inline
getNPSKCarrierFrequencyVector(const std::vector<std::complex<float>> &input,
                              int n = 2)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getNPSKCarrierFrequencyVector()" << std::endl;
#endif

    return getPSD(input, n);
}



/// @brief Vektor mit sich selber konjugiert-komplex multiplizieren (z.B.: AKF).
static std::vector<std::complex<float>> inline
getConj(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getConj(1)" << std::endl;
#endif
    return getConj(input, input);;
}
/// @brief Zwei Vektoren konjugiert-komplex multiplizieren. (z.B.: KKF).
static inline std::vector<std::complex<float>>
getConj(const std::vector<std::complex<float>> &input,
        const std::vector<std::complex<float>> &input2)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getConj(2)" << std::endl;
#endif
    if(input.size() != input2.size() || input.empty())
    {
        throw std::runtime_error("FEHLER libdsss::getConj(2): "
                                 "input.size() != input2.size() || empty");
    }

    std::vector<std::complex<float>> output(input.size());
    std::transform(input.begin(), input.end(),
                   input2.begin(),
                   output.begin(),
                   [](std::complex<float> value, std::complex<float> value2)
                   {return value * std::conj(value2);});

    return output;
}
/// @brief Zwei Vektoren konjugiert-komplex multiplizieren. (z.B.: KKF).
static inline std::vector<std::complex<float>>
getConj(const std::vector<std::complex<float>> &input,
        const std::complex<float> value)
{
    std::vector<std::complex<float>> output(input.size());
    std::complex<float> value_conj = std::conj(value);
    std::transform(input.begin(), input.end(),
                   output.begin(),
                   [value_conj](std::complex<float> sample)
                   {return sample * value_conj;});

    return output;
}



/// @brief Bestimmt die Traegerfrequenz durch potenziern des Eingagnssignals
/// @param input Eingangssignal
/// @return Position der staerksten Frequenz im Vektor
double
getNPSKCarrierFrequency(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getNPSKCarrierFrequency()" << std::endl;
#endif
    std::vector<double> tmp2(getNPSKCarrierFrequencyVector(input));

    std::vector<double>::iterator iter = std::max_element(tmp2.begin(),
                                                          tmp2.end());


    return input.size() - std::distance(tmp2.begin(), iter);
}


/// @brief Zieht die angegebene Frequenz vom Eingangssignal ab
/// @param input Eingangssignal
/// @param frequency_shift Frequenz
/// @return Eingangssignal abzueglich der angegebenen Frequenz
std::vector<std::complex<float>>
shiftFrequency(const std::vector<std::complex<float>> &input,
               int frequency_shift)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::shiftFrequency()" << std::endl;
#endif

    std::vector<std::complex<float>> output(m_fft2->fft(input));

    if(frequency_shift < 0)
    {
        /* Nach rechts rotieren. */
        std::rotate(output.rbegin(),
                    output.rbegin() + std::abs(frequency_shift),
                    output.rend());

    }
    else
    {
        /* Nach rechts rotieren. */
        std::rotate(output.begin(),
                    output.begin() + frequency_shift,
                    output.end());
    }

    return m_fft2->fft(output, 1);
}



/// @brief  Berechnet das Zyklo-Spektrum aus input und addiert es auf base.
/// @param  base Zyklo-Spektrum, bei erster Ausfuehrung mit Null initialisiert
/// @param  input aufzuaddierendes Zyklo-Spektrum
/// @param  exp FFT-Exponent (es muss gelten 2^exp == input.size == base.size)
static inline void
getCylcoSpectrum(std::vector<std::vector<std::complex<float>>> &base,
                 const std::vector<std::complex<float>> &input,
                 uint64_t exp = 6)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getCylcoSpectrum()" << std::endl;
#endif

    FFT2 fft2(exp);
    FFTWindow fftwin;

    getCylcoSpectrum(base, input, fft2, fftwin);
}
static inline void
getCylcoSpectrum(std::vector<std::vector<std::complex<float>>> &base,
                 const std::vector<std::complex<float>> &input,
                 FFT2 &fft2, FFTWindow &fftwin)
{
    std::vector<std::complex<float>> input_fft(
                                fft2.fft(fftwin.hammingWindow(
                                        getNPow<std::complex<float>>(input))));



    // Phasenapassung
//    const float fact = 2.0 * M_PI * i++ / 4;
//    for (size_t j = 0; j < input_fft.size(); j++)
//    {
//        input_fft[j] *= std::exp( std::complex<float>(0, j * fact));
//    }

    std::vector<std::vector<std::complex<float>>> output(input.size());
    for(uint64_t w = 0; w < input_fft.size(); ++w)
    {
        output[w] = Dsss::getConj(input_fft, input_fft[w]);
    }

    if(base.size() != input.size()){base = output;}
    else{Dsss::addToVector(base, output);}
}
/// @brief Eigener Funktionsname, da keine ueberladenen Funktion fuer
///        Threads angenommen werden.
static inline void
getCylcoSpectrumThread(std::vector<std::vector<std::complex<float>>> &base,
                       const std::vector<std::complex<float>> &input,
                       FFT2 &fft2, FFTWindow &fftwin)
{
    getCylcoSpectrum(base, input, fft2, fftwin);
}



/// @brief  Berechnet die zyklisch-verschobene Autokorrelation.
/// @param  Abtastwerte
/// @param  FFT Exponent
/// @return quadratischer Vektor.
std::vector<std::vector<double>>
getCylcoAutoCorr(const std::vector<std::complex<float>> &input,
                 uint64_t exp = 6)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getCylcoAutoCorr()" << std::endl;
#endif
    if(input.empty())
    {
        throw std::runtime_error("FEHLER libdsss::getCylcoAutoCorr(): "
                                 "input.empty()");
    }

    FFT2 fft2(exp);
    std::vector<std::complex<float>> tmp(fft2.fft(input));

    std::vector<std::vector<double>> output(input.size());
    for(uint64_t w = 0; w < tmp.size(); ++w)
    {
        std::rotate(tmp.begin(), tmp.begin() + 1,
                    tmp.end());

        output[w] = getCenteredVector(getLog10(fft2.fft(getConj(tmp), 1)));
    }


    return output;
}



/// @brief Rotiert einen Vektor um die halbe Laenge, so das Null mittig ist.
/// @param input Eingangsvektor (gerade Laenge)
/// @return Zentrierter Eingangsvektor
template <typename T = double>
static std::vector<T>
getCenteredVector(const std::vector<T> &input)
{
    if(input.size() % 2)
    {
        throw std::runtime_error("FEHLER libdsss::getCenteredVector(1): "
                                 "input.size() ist ungerade");
    }

    std::vector<T> output(input.size());
    std::copy(input.begin(), input.begin() + input.size() / 2,
              output.begin() + output.size() / 2);
    std::copy(input.begin() + input.size() / 2, input.end(),
              output.begin());

    return output;
}
/// @brief Rotiert einen Vektor um die halbe Laenge, so das Null mittig ist.
/// @param input Eingangsvektor (gerade Laenge)
/// @return Zentrierter Eingangsvektor
template <typename T = double>
static std::vector<std::vector<T>>
getCenteredVector(const std::vector<std::vector<T>> &input)
{
    if(input.size() % 2)
    {
        throw std::runtime_error("FEHLER libdsss::getCenteredVector(2): "
                                 "input.size() ist ungerade");
    }

    std::vector<std::vector<T>> output(input.size());
    for(uint64_t w = 0; w < input.size() / 2; ++w)
    {
        output[w] = getCenteredVector(input[input.size() / 2 + w]);
        output[input.size() / 2 + w] = getCenteredVector(input[w]);
    }

    return output;
}



/// @brief Gibt das gefensterte Frequenzspektrum aus.
/// @param Eingangsvektor
/// @param
std::vector<double>
getSpectrum(const std::vector<std::complex<float>> &input, int window_index,
            bool log10 = true)
{
    if(input.size() != m_fft_leng)
    {
        throw std::runtime_error("FEHLER libdsss::getSpectrum() != leng");
    }

    if(log10)
    {
        return getLog10(m_fft2->fft(getWindowedSignal(input,
                                                                window_index)));
    }
    return getAbs(m_fft2->fft(getWindowedSignal(input,
                                                                window_index)));
}



/// @brief Resampelt ein Rechtecksignal, es erfolg somit keine Filterung
/// @param input Rechtecksignal
/// @param ratio Samplerate ist / soll
std::vector<int>
resampleSquare(const std::vector<int> &input, double ratio)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::resampleSquare()" << std::endl;
#endif

    double end_of_input = static_cast<double>(input.size() - 1),
           next_sample_pos = 0.0;

    std::vector<int> output;
    output.reserve(input.size());

    while(next_sample_pos < end_of_input)
    {
        double index = 0,
               frac = std::modf(next_sample_pos, &index);

        frac < 0.5 ? (output.push_back(input[static_cast<int>(index)]))
                   : (output.push_back(input[static_cast<int>(index) + 1])) ;

        next_sample_pos += ratio;
    }

    return output;
}



};


class FastXCorr
{
    std::vector<std::complex<float>> _sequence;
    int m_exp;
    FFT2 *_fft2;
    Dsss *_Dsss;

    FastXCorr(int exp)
    {
        m_exp = exp;
        _fft2 = new FFT2(m_exp);
        _Dsss = new Dsss(m_exp);
    }

    ~FastXCorr()
    {
        delete _fft2;
    }

    void setSequence(const std::vector<int> &input)
    {
        std::vector<std::complex<float>> tmp(1 << m_exp,
                                             std::complex<float>(0.0 , 0.0));
        for(uint64_t w = 0; w < input.size(); ++w)
        {
            tmp[w].real(input[w]);
        }

        _sequence = _fft2->fft(tmp);
    }

    std::vector<double> getXCorr(const std::vector<std::complex<float>> &input)
    {
        std::vector<std::complex<float>> tmp(_fft2->fft(input));
        std::transform(tmp.begin(), tmp.end(), _sequence.begin(), tmp.begin(),
                       [](std::complex<float> val_0, std::complex<float> val_1)
                       {return val_0 * std::conj(val_1);});

        return Dsss::getAbs(_fft2->fft(tmp, 1));
    }
};


class CycloSpec
{
    FFT2 *m_fft2;
    FFTWindow m_fftwin;

    uint64_t m_fft_block_leng;

    std::vector<std::vector<std::complex<float>>> m_data;
    std::vector<std::complex<float>> m_buffer;

public:

CycloSpec(int exp)
{
    m_fft_block_leng = 1 << exp;

    m_fft2 = new FFT2(exp);

    m_data.resize(1 << exp);
    for(uint64_t w = 0; w < m_data.size(); ++w)
    {
        m_data[w].resize(1 << exp);
    }
    setDataToZero();
}
~CycloSpec()
{

}


void
setDataToZero(void)
{
    for(auto row : m_data)
    {
        std::fill(row.begin(), row.end(), std::complex<float>(0.0, 0.0));
    }
}


/// @brief Gibt aktuellen
/// @param is_log10 true: 10 * log10() Format (default), false: nur Absolutwerte
/// @return Quadratischer Vektor (2^exp)
std::vector<std::vector<double>>
getResult(bool is_log10 = true)
{
    if(is_log10)
    {
        return Dsss::getLog10(Dsss::getCenteredVector(m_data));
    }
    return Dsss::getAbs(Dsss::getCenteredVector(m_data));
}


/// @brief Fuegt einen Datensatz mit FFT-Laenge hinzu.
/// @param input Laenge muss FFT-Laenge entsprechen.
/// @return 0: hat funktioniert, -1: hat nicht funktioniert
void
addBlock(const std::vector<std::complex<float>> &input)
{
    std::vector<std::complex<float>> input_fft(m_fft2->fft(m_fftwin.hammingWindow(input)));
    std::vector<std::vector<std::complex<float>>> cyclo_tmp(input_fft.size());

    for(uint64_t w = 0; w < cyclo_tmp.size(); ++w)
    {
        cyclo_tmp[w] = Dsss::getConj(input_fft, input_fft[w]);
    }


    return Dsss::addToVector(m_data, cyclo_tmp);
}

};

class GPSCorrelator : Dsss
{

/// Nummer der GPS Satelliten und gleichzeitig
std::vector<std::vector<int>> m_gps_nr =
{
   {1,5}, {2,6}, {3,7}, {4,8}, {0,8}, {1,9}, {0,7}, {1,8}, {2,9}, {1,2}, {2,3},
   {4,5}, {5,6}, {6,7}, {7,8}, {8,9}, {0,3}, {1,4}, {2,5}, {3,6}, {4,7}, {5,8},
   {0,2}, {3,5}, {4,6}, {5,7}, {6,8}, {7,9}, {0,5}, {1,6}, {2,7}, {3,8}
};

/// @return 1023-Bit Satelliten-abhaengigen Gold-Code.
std::vector<int>
getGPSGoldCode(int gps_sat_number)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getGPSGoldCode(1)" << std::endl;
#endif
    return getGPSGoldCode(m_gps_nr[gps_sat_number]);
}

std::vector<int>
getGPSGoldCode(const std::vector<int> &gps_sat_number)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "libdsss::getGPSGoldCode(2)" << std::endl;
#endif
    std::vector<int> result(1023),
                     g1 = {1,1,0,0,1,0,0,0,0,0},
                     g2 = {0,1,0,0,1,0,1,1,0,0};
    for(uint64_t w = 0; w < 1023; ++w)
    {
        g1[9] ^ g2[gps_sat_number[0]] ^ g2[gps_sat_number[1]] ?
                    (result[w] = +1) : (result[w] = -1);


        int fb1 = g1[2] ^ g1[9],
            fb2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9];

        for(uint64_t x = g1.size() - 1; x > 0; --x)
        {
            g1[x] = g1[x -1];
            g2[x] = g2[x -1];
        }

        g1[0] = fb1;
        g2[0] = fb2;
    }

    return result;
}

};


#endif // libdsss_HPP
