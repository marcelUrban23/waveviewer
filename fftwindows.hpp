#ifndef FFTWINDOW_HPP
#define FFTWINDOW_HPP

#include <vector>
#include <complex>

class FFTWindow
{

std::vector<float> lookup_hamming_window,
                   lookup_vonhann_window,
                   lookup_blackman_window,
                   lookup_flattop_window;

const std::vector<std::string> m_windows = {"hamming", "vonhann", "blackman",
                                            "flattop"};

public:

FFTWindow()
{

}

~FFTWindow()
{
    lookup_hamming_window.resize(0);
    lookup_vonhann_window.resize(0);
    lookup_blackman_window.resize(0);
    lookup_flattop_window.resize(0);
}


std::vector<std::string>
getWindows(void) {return m_windows;}


template <typename T = std::complex<float>>
std::vector<T>
getWindowedSignal(const std::vector<T> &input, int index)
{
    return getWindowedSignal(input.data(), input.size(), index);
}
template <typename T = std::complex<float>>
std::vector<T>
getWindowedSignal(const T* input, int leng, int index)
{
    if(index == 0) {return hammingWindow(input, leng);}
    else if(index == 1) {return vonHannWindow(input, leng);}
    else if(index == 2) {return blackmanWindow(input, leng);}
    else {return flattopWindow(input, leng);}
}


template <typename T = std::complex<float>>
std::vector<T>
hammingWindow(const std::vector<std::complex<float>> &input)
{
    return hammingWindow(input.data(), input.size());
}

template <typename T = std::complex<float>>
std::vector<T>
hammingWindow(const T* input, uint64_t leng)
{
    /* Bei der ersten Benutzung initialisieren. */
    if(lookup_hamming_window.size() != leng)
    {
        lookup_hamming_window.resize(leng);

        double constant = M_PI * 2 / (leng - 1);

        for(uint64_t w = 0; w < leng; ++w)
        {
            lookup_hamming_window[w] = static_cast<float>(0.54 - 0.46
                                                         * (cos(constant * w)));
        }
    }

    std::vector<T> output(leng);
    for(uint64_t w = 0; w < leng; ++w)
    {
        output[w] = input[w] * lookup_hamming_window[w];
    }

    return output;
}

template <typename T = std::complex<float>>
std::vector<T>
vonHannWindow(const std::vector<T> &input)
{
    return vonHannWindow(input.data(), input.size());
}
template <typename T = std::complex<float>>
std::vector<T>
vonHannWindow(const T* input, uint64_t leng)
{
    /* Bei der ersten Benutzung initialisieren. */
    if(lookup_vonhann_window.size() != leng)
    {
        lookup_vonhann_window.resize(leng);

        double constant = M_PI * 2 / (leng - 1);

        for(uint64_t w = 0; w < leng; ++w)
        {
            lookup_vonhann_window[w] = static_cast<float>(0.5 - 0.5
                                                         * (cos(constant * w)));
        }
    }

    std::vector<std::complex<float>> output(leng);
    for(uint64_t w = 0; w < leng; ++w)
    {
        output[w] = input[w] * lookup_vonhann_window[w];
    }

    return output;
}

template <typename T = std::complex<float>>
std::vector<T>
blackmanWindow(const std::vector<T> &input)
{
    return vonHannWindow(input.data(), input.size());
}
template <typename T = std::complex<float>>
std::vector<T>
blackmanWindow(const T* input, uint64_t leng)
{
    /* Bei der ersten Benutzung initialisieren. */
    if(lookup_blackman_window.size() != leng)
    {
        double alpha   = 0.16,
               alpha_0 = (1 - alpha) / 2,
               alpha_1 = 1 / 2,
               alpha_2 = alpha / 2;

        lookup_blackman_window.resize(leng);

        double constant = M_PI * 2 / (leng - 1);

        for(uint64_t w = 0; w < leng; ++w)
        {
            lookup_blackman_window[w] = static_cast<float>(
                                           alpha_0
                                           - alpha_1 * (cos(constant * w))
                                           + alpha_2 * (cos(2 * constant * w)));
        }
    }

    std::vector<std::complex<float>> output(leng);
    for(uint64_t w = 0; w < leng; ++w)
    {
        output[w] = input[w] * lookup_blackman_window[w];
    }

    return output;
}


template <typename T = std::complex<float>>
std::vector<T>
flattopWindow(const std::vector<T> &input)
{
    return vonHannWindow(input.data(), input.size());
}
template <typename T = std::complex<float>>
std::vector<T>
flattopWindow(const T* input, uint64_t leng)
{
    /* Bei der ersten Benutzung initialisieren. */
    if(lookup_flattop_window.size() != leng)
    {
        double alpha_0 = 1,
               alpha_1 = 1.93,
               alpha_2 = 1.29,
               alpha_3 = 0.388,
               alpha_4 = 0.028;

        lookup_flattop_window.resize(leng);

        double constant = M_PI * 2 / (leng - 1);

        for(uint64_t w = 0; w < leng; ++w)
        {
            lookup_flattop_window[w] = static_cast<float>(
                                           alpha_0
                                           - alpha_1 * (cos(constant * w))
                                           + alpha_2 * (cos(2 * constant * w))
                                           + alpha_3 * (cos(3 * constant * w))
                                           + alpha_4 * (cos(4 * constant * w)));
        }
    }

    std::vector<std::complex<float>> output(leng);
    for(uint64_t w = 0; w < leng; ++w)
    {
        output[w] =  input[w] * lookup_flattop_window[w];
    }

    return output;
}

};

#endif // FFTWINDOW_HPP
