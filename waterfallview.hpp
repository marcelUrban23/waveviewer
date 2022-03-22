#ifndef WATERFALLVIEW_HPP
#define WATERFALLVIEW_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <QMainWindow>
#include <QTextEdit>
#include <QLineEdit>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QComboBox>
#include <QGroupBox>
#include <QScrollBar>
#include <QKeyEvent>
#include <QPushButton>
#include <QCheckBox>
#include <QFile>
#include <QFileDialog>
#include <QLabel>
#include <QPainter>
#include <QPixmap>
#include <QPen>
#include <QMessageBox>
#include <QPaintEvent>

//#define DEBUG_FUNCTION_CALL

/// @brief Zeichnet ein Graues Wasserfalldiagramm
class WaterfallView : public QWidget
{
Q_OBJECT

std::vector<std::vector<float>> m_data;
std::vector<QPen> m_qpens;
QPixmap  m_qpm_spec;

bool m_is_horizontal = true;

public:


/// @brief Erstellt ein Wasserfalldiagramm.
/// @param title Fenstertitel
/// @param fft_exp Aufloesung der FFT, nachtraeglich nicht veraenderbar.
WaterfallView(QString title = "Hier_koennte_Ihr_Name_stehen",
              QWidget *parent = nullptr):
              m_is_horizontal(true)
{
    std::cerr << "WaterfallView::WaterfallView()" << std::endl;

    setWindowTitle(title);



    // Stifte mit Grauwerten
    for(int w = 0; w < 255; ++w)
    {
        m_qpens.push_back(QPen(QColor(0, 0, 0, w)));
    }
}

~WaterfallView()
{
    m_data.resize(0);
    m_qpens.resize(0);
}


/// @brief Schiebt Abtastwerte in die gepufferte Verarbeitungskette.
/// @param input Vektor wird je nach Laenge in FFT-Bloecke aufgeteilt und
///        dem Spektrogramm angefuegt.
/// @param is_fft_transformed gibt an, ob der Eingangsvektor bereits fft ist.
void
pushSamples(const std::vector<std::vector<std::complex<float>>> &input)
{
#ifndef DEBUG_FUNCTION_CALL
    std::cerr << "WaterfallView::pushSamples(1)" << std::endl;
#endif

    m_data.resize(input.size());
    for(uint64_t w = 0; w < m_data.size(); ++w)
    {
        m_data[w].resize(input[w].size());
        std::transform(input[w].begin(), input[w].end(), m_data[w].begin(),
                       [](std::complex<float> sample) -> float
                       {return static_cast<float>(std::abs(sample));});
    }
}
void
pushSamples(const std::vector<std::vector<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaterfallView::pushSamples(2)" << std::endl;
#endif
    m_data = input;

    update();
}


private:

///// @brief Reagiert auf Groessenaenderung
//void
//resizeEvent(QResizeEvent *qre)
//{
//    std::cerr << "WaterfallView::resizeEvent()" << std::endl;

//    m_qpm_spec = m_qpm_spec.scaled(this->size().width(), this->size().height());
//    QPainter qp_main(this);
//    qp_main.drawPixmap(0, 0, m_qpm_spec);
//}


/// @brief Zeichnet
void
paintEvent(QPaintEvent *qpe)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaterfallView::paintEvent()" << std::endl;
#endif
    if(m_data.empty()){return;}

    float subst_w = static_cast<float>(this->size().width())
                  / static_cast<float>(m_data.size());
    float subst_x = static_cast<float>(this->size().height())
                  / static_cast<float>(m_data[0].size());

    std::cerr << "subst_w: " << subst_w << std::endl;
    std::cerr << "subst_x: " << subst_x << std::endl;

    std::vector<std::vector<float>> tmp_data(this->size().width());
    for(uint64_t w = 0; w < tmp_data.size(); ++w)
    {
        tmp_data[w] = std::vector<float>(this->size().height(), 0.0);
    }
    if(subst_w < 1.0) // Horizontal Interpolieren
    {
        for(uint64_t w = 0; w < m_data.size(); ++w)
        {
            if(subst_x < 1.0) // Vertikal Interpolieren
            {
                for(uint64_t x = 0; x < m_data[0].size(); ++x)
                {
                    tmp_data[static_cast<uint64_t>(w * subst_w)]
                            [static_cast<uint64_t>(x * subst_x)]
                            += m_data[w][x];
                }
            }
            else // Vertikal Expandieren
            {
                for(uint64_t x = 0;
                    static_cast<uint64_t>(x / subst_x) < m_data[0].size(); ++x)
                {
                    tmp_data[static_cast<uint64_t>(w * subst_w)]
                            [x]
                    += m_data[w]
                             [static_cast<uint64_t>(x / subst_x)];
                }
            }
        }
    }
    else // Horizontal expandieren
    {
        for(uint64_t w = 0;
            static_cast<uint64_t>(w / subst_w) < m_data.size(); ++w)
        {
            if(subst_x < 1.0) // Vertikal Inter2polieren
            {
                for(uint64_t x = 0; x < m_data[0].size(); ++x)
                {
                    tmp_data[w]
                            [static_cast<uint64_t>(x * subst_x)]
                   += m_data[static_cast<uint64_t>(w / subst_w)]
                            [x];
                }
            }
            else // Vertikal Expandieren
            {
                for(uint64_t x = 0;
                    static_cast<uint64_t>(x / subst_x) < m_data[0].size(); ++x)
                {
                    tmp_data[w][x] +=m_data[static_cast<uint64_t>(w / subst_w)]
                                           [static_cast<uint64_t>(x / subst_x)];
                }
            }
        }
    }

    float amp_max = tmp_data[0][0],
          amp_min = amp_max;

    for(auto row : tmp_data)
    {
        float max = *std::max_element(row.begin(), row.end());
        float min = *std::min_element(row.begin(), row.end());

        if(max > amp_max){amp_max = max;}
        if(min < amp_min){amp_min = min;}
    }
    // Ermittelt die Schrittweite der Grauwerte in Abhaengigkeit der Dynamik.
    float grey_stepwidth = (amp_max - amp_min) /  254.0;

    std::cerr << "amp_max: " << amp_max << ", amp_min: " << amp_min
              << ",grey_stepwidth: " << grey_stepwidth << std::endl;

    // Wasserfallbild vertialer Verlauf
    if( ! m_is_horizontal)
    {

    }
    m_qpm_spec = QPixmap(tmp_data.size(), tmp_data[0].size());
    m_qpm_spec.fill(Qt::white);
    QPainter qp_spec(&m_qpm_spec);



    for(uint64_t w = 0; w < tmp_data.size(); ++w)
    {
        for(uint64_t x = 0; x < tmp_data[w].size(); ++x)
        {
            // Farbe entsprechend des Betrages waehlen
            qp_spec.setPen(m_qpens[static_cast<uint32_t>((  tmp_data[w][x]
                                                          - amp_min)
                                                          / grey_stepwidth)]);

            // Frequenz auftragen
            qp_spec.drawPoint(w, x);
        }
    }

    // Beispiel dazu ansehen, um Fehler zu erkennen
//    m_qpm_spec = m_qpm_spec.scaled(this->size().width(),
//                                   this->size().height());

    QPainter qp_main(this);
    qp_main.drawPixmap(0, 0, m_qpm_spec);
}

};
#endif // WATERFALLVIEW_HPP
