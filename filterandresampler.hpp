#ifndef FILTERANDRESAMPLER_HPP
#define FILTERANDRESAMPLER_HPP

#include <vector>
#include <complex>

#include "libdsss.hpp"

class FilterAndResampler
{

/// @brief Berechnet ein Filterfenster wobei der erste Teil der
///        Durchlassbereich, der zweite Teile der Uebergangsberich und der
///        dritte Teil der Daempfungsbereich ist.
/// @param leng Laenge der zu filternden Signalbloecke
/// @param passband Durchlassbereich vor dem Uebergangsbereich
/// @param transition Uebergangsbereich vor dem Stopband (Steilheit des Filters)
/// @param stopband_attenuation Daempfung des Stopbandes
/// @return Filterbins fuer den Frequenzbereich
template <typename T = std::complex<float>>
std::vector<T>
getLowPassWindow(uint64_t leng, double passband, double transition,
                 double stopband_attenuation = -60.0)
{
    transition = std::min(std::max(1e-6, std::abs(transition)), 0.99);
    if(leng <= 1)
    {
        throw std::invalid_argument("libdsss::getLowPassWindow: "
                                    "eng <= 1");
    }
    if(stopband_attenuation >= 1.0)
    {
        throw std::invalid_argument("libdsss::getLowPassWindow: "
                                    "stopband_attenuation >= 1.0");
    }
    double pass_factor = 1.0,
           stop_factor = std::pow(10.0, stopband_attenuation / 10.0);

    // Das Fenster wird fuer ein Seitenband berechnet und spaeter gespiegelt
    std::vector<double> window(leng / 2, pass_factor);

    uint64_t i = std::max<uint64_t>(0, static_cast<int>(std::ceil(passband
                                                        * window.size())));
    for(; i < window.size(); ++i)
    {
        double x = static_cast<double>(i) / (window.size() - 1) - passband;
        if(x > transition) // Stopband
        {
            window[i] = stop_factor;
        }
        else // Uebergangsbereich (Steilheit des Filters)
        {
            window[i] = (0.5 + 0.5 * std::cos(M_PI * x / transition))
                        * (pass_factor - stop_factor) + stop_factor;
        }
    }

    std::vector<T> window_T;
    window_T.reserve(leng);

    std::transform(window.begin(), window.end(),
                   std::back_inserter(window_T),
                   [](double value){return static_cast<T>(value);});

    if(leng % 2) window_T.push_back(window.back());

    std::transform(window.rbegin(), window.rend(),
                   std::back_inserter(window_T),
                   [](double value){return static_cast<T>(value);});

    return window_T;
}


/// @brief Interpoliert einen Punkt zwischen
/// @param position relative Position zwischen zwei Punkten
/// @return interpolierter Punkt zwischen zwei Punkten bzw. wenn
///         kein Vorgaenger oder Nachfolger: letzter oder erster
template <typename T = std::complex<float>>
T
interpolateLinear(const std::vector<T> &input, double position)
{
    if(position <= 0.0) return input[0];
    if(position > input.size() - 1) return input.back();

    double ipos,
           fpos = std::modf(position, &ipos);

    uint64_t i = static_cast<uint64_t>(ipos);

    return    input[i] * static_cast<T>(1.0 - fpos)
            + input[i + 1] * static_cast<T>(fpos);
}

/// @brief
template <typename T = std::complex<float>>
T
interpolCubic(const std::vector<T> &input, double pos)
{
    if(input.size() < 3 || pos < 1.0 || pos >= input.size() - 2)
    {
        return interpolateLinear<T> (input, pos);
    }

    double ipos,
           fpos = std::modf(pos, &ipos);

    uint64_t i = static_cast<uint64_t>(ipos);

    T a_0 = -static_cast<T>(0.5) * (input[i - 1] - input[i + 2])
            + static_cast<T>(1.5) * (input[i] - input[i + 1]),
      a_1 = input[i - 1] - static_cast<T>(2.5) * (input[i] - input[i + 1])
            - static_cast<T>(0.5) * input[i + 2],
      a_2 = -static_cast<T>(0.5) * (input[i - 1] - input[i + 1]),
      a_3 = input[i];

    return   a_0 * static_cast<T>(fpos * fpos * fpos)
           + a_1 * static_cast<T>(fpos * fpos)
           + a_2 * static_cast<T>(fpos)
           + a_3;
}


public:

/// @brief Verschiebt das Signal um shift
/// @param samplerate Abtastfrequenz des Eingagnssignals
/// @param shift Frequenz, um die das Signal verschoben werden soll
/// @param verschobenes Eingangssignal
template <typename T = std::complex<float>>
std::vector<T>
shiftFrequency(const std::vector<T> &input, double samplerate, double shift)
{
    return shiftFrequency(input, shift / samplerate);
}
template <typename T = std::complex<float>>
std::vector<T>
shiftFrequency(const std::vector<T> &input, double shift_ratio)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "FilterAndResampler::shiftFrequency()" << std::endl;
#endif

    FFT2 fft(Dsss::getNextPow(input.size()));

    std::vector<T> data = input;
    data.resize(fft.getLeng(), T(0.0, 0.0));

    std::vector<std::complex<float>> output(fft.fft(data));

    int64_t shift = static_cast<int64_t>(shift_ratio
                                         * static_cast<double>(fft.getLeng()));

    std::cerr << "shift_ratio: " << shift_ratio << ", fft->leng: "
              << fft.getLeng() << std::endl;

    if(shift < 0)
    {
        /* Nach rechts rotieren. */
        std::rotate(output.rbegin(),
                    output.rbegin() + std::abs(shift),
                    output.rend());
    }
    else
    {
        /* Nach rechts rotieren. */
        std::rotate(output.begin(),
                    output.begin() + shift,
                    output.end());
    }

    output = fft.fft(output, 1);
    output.resize(input.size());

    return output;
}

public:

FilterAndResampler(void)
{
}


template <typename T = std::complex<float>>
std::vector<T>
getFilteredSignal(const std::vector<T> &input, double passband,
                                               double transisiton_width)
{
    FFT2 fft(Dsss::getNextPow(input.size()));



    std::vector<std::complex<float>> lowpasswin(getLowPassWindow(fft.getLeng(),
                                                passband, transisiton_width));

    std::vector<T> data = input;
    data.resize(fft.getLeng(), T(0.0, 0.0));

    std::vector<std::complex<float>> input_fft(fft.fft(data));

    std::transform(input_fft.begin(), input_fft.end(), lowpasswin.begin(),
                   input_fft.begin(),
                   std::multiplies<std::complex<float>>());

    std::vector<std::complex<float>> output(fft.fft(input_fft, 1));

    // Amplitude anpassen
//    std::vector<double> cum_sum = getCumulativeSum(input);

//    for(uint64_t w = 0; w < input.size(); ++w)
//    {
//        output *= static_cast<T>(SIZE / cum_sum[w]);
//    }

    output.resize(input.size());

    return output;
}

/// @brief Resampelt den Eingangsvektor entsprechend der Rate
/// @param input Vektor mit Abtastwerten
/// @param ratio Abtastrate: Eingang / Ausgang
/// @return resampelter Vektor mit ungleicher Groesse bez. des Eingangsvektors
template <typename T = std::complex<float>>
std::vector<T>
resampleCubic(const std::vector<T> &input, double ratio)
{
    double start = 0.0;

    std::vector<T> output;
    output.reserve(static_cast<uint64_t>(input.size() / ratio) + 1);

    while(start < input.size() && start >= 0.0)
    {
        output.push_back(interpolCubic(input, start));
        start += ratio;
    }

    const T norm = static_cast<T>(Dsss::getVariance(input)
                                  / Dsss::getVariance(output));
    std::transform(output.begin(), output.end(),output.begin(),
                   [norm](T &sample)
                   {return sample * norm;});


    return output;
}

};

#endif // FILTERANDRESAMPLER_HPP
