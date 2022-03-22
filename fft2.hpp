#ifndef FFT2_HPP
#define FFT2_HPP

#include <stdint.h>
#include <stdlib.h>
#include <iostream> // std::cerr

#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>

//#define DEBUG

class FFT2
{
    std::vector<std::complex<float>> m_omaga_sorted;

    std::vector<uint64_t> m_result_index_sorted;

    uint64_t m_fft2_leng,
             m_exp;

public:


/* Konstruktor */
FFT2(uint64_t exponent)
{
#ifdef DEBUG
    fprintf(stderr, "INFO fft2::fft2()\n");
#endif

    double   rad  = 0.0;

    uint64_t w    = 0;
    m_exp         = exponent;
    m_fft2_leng   = 1 << m_exp;


    /* Speicher reservieren. */
    m_omaga_sorted.resize(m_fft2_leng);

    m_result_index_sorted.resize(m_fft2_leng);

    std::vector<std::complex<float>> omega(m_fft2_leng);



    /* Quant berechnen. */
    rad = (2.0 * M_PI) / m_fft2_leng;

    /* Erstellen der Drehfaktor (twiddle factors) -Liste. */
    for(w = 0; w < m_fft2_leng; ++w)
    {
        omega[w].real(  cos(rad * w));
        omega[w].imag( -sin(rad * w));
    }

    /* Erstellen der Zuordnungsliste fuer die Ergebnisse. */
    for(w = 0; w < m_fft2_leng; ++w)
    {
        m_result_index_sorted[w] = reverseBits(w, m_exp);
    }

    /* Erstellen der sortierten Drehfaktor-Index-Liste. */
    for(w = 0; w < (m_fft2_leng / 2); ++w)
    {
        m_omaga_sorted[w] = omega[m_result_index_sorted[2 * w]];
    }


    return;
}


/* Destruktor */
~FFT2()
{

    return;
}

uint64_t
getLeng(void)
{
    return m_fft2_leng;
}


/// @brief Diese Funktion berechnet die Fouriertransformierte der
///        Eingangswerte in Form eines Vektors mit der gleichen
///        Anzahl Ausgangswerte.
/// @param input
/// @param is_inverse = 0: Hintransformation, 1: Ruecktransformation
/// @return Fouriertransformierte
std::vector<std::complex<float>>
fft(std::vector<std::complex<float>> input, int64_t is_inverse = 0)
{
    return fft(input.data(), input.size());
}
std::vector<std::complex<float>>
fft(std::complex<float> *input, uint64_t leng, int64_t is_inverse = 0)
{
    if(leng != m_fft2_leng)
    {
        throw std::runtime_error("FEHLER FFT2::fft: input.size() != m_fft2_leng");
    }
    if(is_inverse){return fft_backward(input, leng);}
    else          {return fft_forward(input, leng);}
}

private:

std::vector<std::complex<float>>
fft_backward(std::complex<float> *input, uint64_t leng)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "FFT2::fft_backward()" << std::endl;
#endif

    /* Speicher fuer die Rueckgabewerte reservieren. */
    std::vector<std::complex<float>> output(m_fft2_leng),
                                     &omega_sorted = m_omaga_sorted;

    std::vector<uint64_t> &result_index_sorted = m_result_index_sorted;

    std::complex<float> second_note;
    float real, imag;
    uint64_t phase_index,
             interval,
             blocks,
             offset,
             start,
             n, l, k;

    std::complex<float> frq = omega_sorted[0];


    for(l = 0; l < m_exp; ++l)
    {
        interval = ((m_fft2_leng / 2) >> l);

        blocks = (1 << l);

        phase_index = 0;

        for(k = 0; k < blocks; ++k)
        {
            offset = k * (interval << 1);

            frq = omega_sorted[phase_index];

            for(n = 0; n < interval; ++n)
            {
                start = offset + n;
                second_note = input[start + interval];

                real =  second_note.imag() *  frq.imag()
                     +  second_note.real() *  frq.real();

                imag =  second_note.real() * -frq.imag()
                     +  second_note.imag() *  frq.real();

                input[start + interval].real(  input[start].real()
                                              - real);

                input[start + interval].imag(  input[start].imag()
                                              - imag);

                input[start].real(input[start].real() + real);
                input[start].imag(input[start].imag() + imag);
            }

            ++phase_index;
        }
    }

    /* Ausgangswerte sortieren und normieren. */
    for(k = 0; k < output.size(); ++k)
    {
        output[k] = input[result_index_sorted[k]];
        output[k] /= m_fft2_leng;
    }


    return output;
}


std::vector<std::complex<float>>
fft_forward(std::complex<float> *input, uint64_t leng)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "FFT2::fft_forward()" << std::endl;
#endif

    /* Speicher fuer die Rueckgabewerte reservieren. */
    std::vector<std::complex<float>> output(m_fft2_leng),
                                     &omega_sorted = m_omaga_sorted;

    std::vector<uint64_t> &result_index_sorted = m_result_index_sorted;

    std::complex<float> second_note;
    float real, imag;
    uint64_t phase_index,
             interval,
             blocks,
             offset,
             start,
             n, l, k;

    std::complex<float> frq = omega_sorted[0];


    for(l = 0; l < m_exp; ++l)
    {
        interval = ((m_fft2_leng / 2) >> l);

        blocks = (1 << l);

        phase_index = 0;

        for(k = 0; k < blocks; ++k)
        {
            offset = k * (interval << 1);

            frq = omega_sorted[phase_index];

            for(n = 0; n < interval; ++n)
            {
                start = offset + n;
                second_note = input[start + interval];

                real = -second_note.imag() * frq.imag()
                     +  second_note.real() * frq.real();

                imag =  second_note.real() * frq.imag()
                     +  second_note.imag() * frq.real();


                input[start + interval].real(  input[start].real()
                                              - real);

                input[start + interval].imag(  input[start].imag()
                                              - imag);

                input[start].real(input[start].real() + real);
                input[start].imag(input[start].imag() + imag);
            }

            ++phase_index;
        }
    }

    /* Ausgangswerte sortieren. */
    for(k = 0; k < output.size(); ++k)
    {
        output[k] = input[result_index_sorted[k]];
    }


    return output;
}

/// @brief dreht Bitreihenfolge um.
/// @param bits Bits, deren Reihenfolge umgedreht bzw. gespiegelt werden soll
/// @param bits_leng doppelte Laenge der Spiegelachse
/// @return umgekehrte Bitreihenfolge
uint64_t
reverseBits(uint64_t bits, uint64_t bits_leng)
{
    uint64_t bits_reversed = 0; /* nimmt Bits in umgedrehter Folge auf */

    for(uint64_t w = 0; w < bits_leng; ++w)
    {
        /* Je Schleifendurchgang wird die Stelle "w" im Bitmuster an Position */
        /* des LSB (ganz links) geschoben und verUNDet. Danach wird diese */
        /* entsprechend nach links geschoben. */

        /* Beispliel: "bits": 110101 (bits_leng: 6): */
        /* erster   Durchgang: (((110101 >> (6 - 1 - 0) & 1) << 0 = 000001 */
        /* zweiter  Durchgang: (((110101 >> (6 - 1 - 1) & 1) << 1 = 000011 */
        /* dritter  Durchgang: (((110101 >> (6 - 1 - 2) & 1) << 2 = 000011 */
        /* vierter  Durchgang: (((110101 >> (6 - 1 - 3) & 1) << 3 = 001011 */
        /* fuenf    Durchgang: (((110101 >> (6 - 1 - 4) & 1) << 4 = 001011 */
        /* sechster Durchgang: (((110101 >> (6 - 1 - 5) & 1) << 5 = 101011 */

        bits_reversed |= ((bits >> (bits_leng - 1 - w)) & 1) << w;
    }


    return bits_reversed;
}

};
#endif //FFT2_HPP
