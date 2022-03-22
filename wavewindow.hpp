#ifndef WAVEWINDOW_HPP
#define WAVEWINDOW_HPP

#include <QMainWindow>
#include <QLineEdit>

#include <QMessageBox>
#include <iostream>
#include <numeric>

//#include <QtDataVisualization>

#include "libdsss.hpp"
#include "wavedata.hpp"
#include "waterfallview.hpp"


//#include "extlib/QPlot3D.h"
#include "wavepainter.hpp"



/// @brief Dient primaer der Aufbereitung und Anzeige der Daten
///        Gebaut, um in WaveViewer eingebetet werden.
///        Marker zum aussmessen von Absolutwerten und Differenzen
///        Interpretation der Abtastwerte hinsichtlich verschiedener Darstellung
///        z.B. Real, Imaginaer, Phase, Spektrum, Leistungsdichtespektrum (x^n)
///             Autokorrelation
///
class WaveWindow : public QWidget
{
Q_OBJECT

WaterfallView *m_waterfall;

FFTWindow *m_fftwin;

WavePainter *m_plotter;
QLabel      *m_ql_amount_of_samples,
            *m_ql_samp_rate;
QComboBox   *m_qle_fft_exp,
            *m_qcb_views,
            *m_qcb_windows;
QCheckBox   *m_qcb_window_disabled;
QVBoxLayout *m_qvbl_main;
QComboBox   *m_qcb_waterfall;

QLineEdit   *m_qle_marker_one,
            *m_qle_marker_two,
            *m_qle_marker_diff;

int m_view_index = 0;

std::vector<std::complex<float>> m_data;
const std::vector<std::string> m_views {"Real",
                                        "Imag",
                                        "Phase",
                                        "Spektrum",
                                        "AM-Detekt",
                                        "Leistungsdichtespektrum",
                                        "Leistungsdichtespektrum(x(t)^2",
                                        "Leistungsdichtespektrum(x(t)^4)",
                                        "Autokorrelation",
                                        "Wasserfall"};

bool m_was_waterfall,
     m_is_centered,
     m_is_log10,
     m_is_zoom_to_marker;

uint64_t m_x_offset;

public:

WaveWindow(void) :
    m_was_waterfall(false),
    m_is_centered(true),
    m_is_log10(true),
    m_is_zoom_to_marker(false),
    m_x_offset(0)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::WaveWindow()" << std::endl;
#endif
    createGui();

    setContextMenuPolicy(Qt::CustomContextMenu);

}
~WaveWindow(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::~WaveWindow()" << std::endl;
#endif

    m_data.resize(0);
    delete m_waterfall;
    delete m_fftwin;
    delete m_plotter;
    delete m_ql_amount_of_samples;
    delete m_ql_samp_rate;
    delete m_qle_fft_exp;
    delete m_qcb_views;
    delete m_qcb_windows;
    delete m_qcb_window_disabled;
    delete m_qvbl_main;
    delete m_qcb_waterfall;
}


private:


void createGui(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::createGui()" << std::endl;
#endif
    m_qcb_window_disabled = new QCheckBox;
    m_qcb_window_disabled->setChecked(false);

    QPushButton *qpb_zoom = new QPushButton("zoom");
    connect(qpb_zoom, SIGNAL(clicked(bool)),
            this       , SLOT(zoom()));

    QPushButton *qpb_rezoom = new QPushButton("rezoom");
    connect(qpb_rezoom, SIGNAL(clicked(bool)),
            this       , SLOT(rezoom()));

    QPushButton *qpb_find_max = new QPushButton("Maximalwert");
    connect(qpb_find_max, SIGNAL(clicked(bool)),
            this        , SLOT(findMax()));

    m_ql_amount_of_samples = new QLabel("0");

    m_ql_samp_rate = new QLabel("1");

    m_qcb_views = new QComboBox;
    for(auto view : m_views)
    {
        m_qcb_views->addItem(QString::fromStdString(view));
    }
    connect(m_qcb_views, SIGNAL(currentIndexChanged(int)),
            this       , SLOT(plot()));

    // FFT Exponenten eintragen
    m_qle_fft_exp = new QComboBox;
    for(uint64_t w = 2; w < 28; ++w)
    {
        m_qle_fft_exp->addItem(QString::number(w));
    }
    m_qle_fft_exp->setCurrentIndex(11);
    connect(m_qle_fft_exp, SIGNAL(currentIndexChanged(int)),
            this         , SLOT(plot()));

    m_qcb_windows = new QComboBox;
    m_fftwin = new FFTWindow;
    for(auto win : m_fftwin->getWindows())
    {
        m_qcb_windows->addItem(QString::fromStdString(win));
    }
    connect(m_qcb_windows, SIGNAL(currentIndexChanged(int)),
            this         , SLOT(plot()));

    m_qcb_waterfall = new QComboBox;
    m_qcb_waterfall->addItem("Horizontal");
    m_qcb_waterfall->addItem("Vertikal");

    m_waterfall = new WaterfallView();

    m_plotter = new WavePainter;

    QHBoxLayout *qhbl_marker = new QHBoxLayout;
    qhbl_marker->addWidget(new QLabel("m_1: "));
    m_qle_marker_one = new QLineEdit;
    m_qle_marker_one->setMinimumWidth(100);
    connect(m_qle_marker_one, SIGNAL(editingFinished()),
            this            , SLOT(setMarkerOne()));
    connect(m_plotter       , SIGNAL(markerChanged()),
            this            , SLOT(markerWasChanged()));
    qhbl_marker->addWidget(m_qle_marker_one);

    qhbl_marker->addWidget(new QLabel("m_2: "));
    m_qle_marker_two = new QLineEdit;
    m_qle_marker_two->setMinimumWidth(100);
    connect(m_qle_marker_two, SIGNAL(editingFinished()),
            this            , SLOT(setMarkerTwo()));
    connect(m_plotter       , SIGNAL(markerChanged()),
            this            , SLOT(markerWasChanged()));
    qhbl_marker->addWidget(m_qle_marker_two);

    qhbl_marker->addWidget(new QLabel("Diff: "));
    m_qle_marker_diff = new QLineEdit;
    m_qle_marker_diff->setMinimumWidth(100);
    qhbl_marker->addWidget(m_qle_marker_diff);
    connect(m_qle_marker_diff, SIGNAL(editingFinished()),
            this             , SLOT(markerDiffChanged()));


    QHBoxLayout *qhbl_settings = new QHBoxLayout;
    qhbl_settings->setAlignment(Qt::AlignLeft);

    qhbl_settings->addWidget(new QLabel("Fenster einfrieren: "));
    qhbl_settings->addWidget(m_qcb_window_disabled);

    qhbl_settings->addWidget(qpb_zoom);
    qhbl_settings->addWidget(qpb_rezoom);
    qhbl_settings->addWidget(qpb_find_max);

    qhbl_settings->addWidget(new QLabel("Abtastwerte: "));
    qhbl_settings->addWidget(m_ql_amount_of_samples);

    qhbl_settings->addWidget(new QLabel("Samplerate: "));
    qhbl_settings->addWidget(m_ql_samp_rate);

    qhbl_settings->addWidget(new QLabel("Ansicht: "));
    qhbl_settings->addWidget(m_qcb_views);

    qhbl_settings->addWidget(new QLabel("FFT-Exponent: "));
    qhbl_settings->addWidget(m_qle_fft_exp);

    qhbl_settings->addWidget(m_qcb_windows);

    qhbl_settings->addWidget(m_qcb_waterfall);

    qhbl_settings->addLayout(qhbl_marker);

    QPushButton *qpb_close = new QPushButton("schliessen");
    qhbl_settings->addWidget(qpb_close);
    connect(qpb_close, SIGNAL(clicked()),
            this     , SLOT(close()));

    m_qvbl_main = new QVBoxLayout;
    m_qvbl_main->addLayout(qhbl_settings);
    m_qvbl_main->addWidget(m_plotter);
    m_qvbl_main->addWidget(m_waterfall);
    m_waterfall->setHidden(true);

    this->setLayout(m_qvbl_main);
}




void calculateAnddDrawSignal(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::plot()" << std::endl;
#endif
    if(m_data.empty() || m_qcb_window_disabled->isChecked()){return;}
    int view_index = m_qcb_views->currentIndex();

    if(view_index < 9 && m_was_waterfall)
    {
        m_waterfall->setHidden(true);
        m_plotter->setHidden(false);
        m_was_waterfall = false;
    }
    else if(view_index == 9 && ! m_was_waterfall)
    {
        m_plotter->setHidden(true);
        m_waterfall->setHidden(false);
        m_was_waterfall = true;
    }

    std::vector<std::complex<float>> data(m_data);

    // Hier kommen die zu Zeichnenden (y-)Werte
    std::vector<double> magnitudes(data.size());

    // Auswahl der Darstellung
    if(view_index == 0) // Realteil
    {
        std::transform(data.begin(), data.end(), magnitudes.begin(),
                       [](std::complex<float> sample) -> double
                       {return std::real(sample);});
    }
    else if(view_index == 1) // Imaginaerteil
    {
        std::transform(data.begin(), data.end(), magnitudes.begin(),
                       [](std::complex<float> sample) -> double
                       {return std::imag(sample);});
    }
    else if(view_index == 2)// Phase
    {
        std::transform(data.begin(), data.end(), magnitudes.begin(),
                       [](std::complex<float> sample) -> double
                       {return std::arg(sample);});
    }
    else if(view_index < 9)
    {
//        m_plotter->xAxis->setLabel("Frequenz [Hz]");

        // Daten der FFT Groesse anpassen, ggf. Rest mit Nullen auffuellen.
        data.resize(1 << (m_qle_fft_exp->currentIndex() + 2),
                                              std::complex<float>(0.0, 0.0));

        // + 2 um vonm Index auf den richtigen fft-Exponenten
        Dsss dsss(m_qle_fft_exp->currentIndex() + 2);

        // Auswahl der
        if(view_index == 3) // Frequenzspektrum
        {
            magnitudes = dsss.getSpectrum(data, m_qcb_windows->currentIndex(),
                                          false);
        }
        else if(view_index == 4) // AM-Detekt von x(t)
        {
            magnitudes = dsss.getAMDetect(data, m_qcb_windows->currentIndex());
        }
        else if(view_index == 5) // Leistungsdichtespektrum von x(t)
        {
            magnitudes = dsss.getPSD(data, 1, m_qcb_windows->currentIndex());
        }
        else if(view_index == 6) // Leistungsdichtespektrum von x(t)^2
        {
            magnitudes = dsss.getPSD(data, 2, m_qcb_windows->currentIndex());
        }
        else if(view_index == 7) // Leistungsdichtespektrum von x(t)^4
        {
            magnitudes = dsss.getPSD(data, 4, m_qcb_windows->currentIndex());
        }
        else if(view_index == 8) // Autokorrelation
        {
            magnitudes = Dsss::getAbs(dsss.getAutoCorr(data));
        }

        if(m_is_centered) magnitudes = Dsss::getCenteredVector(magnitudes);
        if(m_is_log10) magnitudes = Dsss::getLog10(magnitudes);
    }

    if(m_is_zoom_to_marker)
    {
        uint64_t start = std::min(m_plotter->getMarkerOnePos(),
                                  m_plotter->getMarkerTwoPos());
        uint64_t stop  = std::max(m_plotter->getMarkerOnePos(),
                                  m_plotter->getMarkerTwoPos());

        start = std::min(start, magnitudes.size() - 1);
        stop  = std::min(stop , magnitudes.size() - 1);

        m_plotter->setXOffset(m_x_offset + start);
        m_plotter->setMagnitudes(magnitudes.begin() + start,
                                 magnitudes.begin() + stop);
    }
    else
    {
        m_plotter->setXOffset(m_x_offset + 0);
        m_plotter->setMagnitudes(magnitudes);
    }

    if(view_index > 8) // Wassserfalldiagramm
    {
        FFT2 fft(m_qle_fft_exp->currentIndex() + 2);

        std::vector<std::vector<float>> sklar;
        for(uint64_t w = 0; w < m_data.size() / fft.getLeng(); ++w)
        {
            sklar.push_back(Dsss::getCenteredVector<float>(Dsss::getLog10<float>(
                                    fft.fft(
                                        m_fftwin->getWindowedSignal(
                                              &(data.data())[w * fft.getLeng()],
                                                             fft.getLeng(),
                                              m_qcb_windows->currentIndex()),
                                                         fft.getLeng()))));
        }

        m_waterfall->pushSamples(sklar);
    }

    update();
}

void close(void)
{
    this->~WaveWindow();
}


private slots:

void plot(void)
{
    calculateAnddDrawSignal();
}

void zoom(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::zoom()" << std::endl;
#endif
    m_is_zoom_to_marker = true;
    plot();
}
void rezoom(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::rezoom()" << std::endl;
#endif
    m_is_zoom_to_marker = false;
    plot();
}


void
setMarkerOne(void)
{
    m_plotter->setMarkerOne(m_qle_marker_one->text());
    m_qle_marker_diff->setText(QString::number(
                                   m_qle_marker_two->text().toULongLong() -
                                   m_qle_marker_one->text().toULongLong()));
}


void
setMarkerTwo(void)
{
    m_plotter->setMarkerTwo(m_qle_marker_two->text());
    m_qle_marker_diff->setText(QString::number(
                                   m_qle_marker_two->text().toULongLong() -
                                   m_qle_marker_one->text().toULongLong()));
}



/// @brief Aendert sich die Differenz, so aendert sich auch nur Marker
void
markerDiffChanged(void)
{
    m_plotter->setMarkerTwo(QString::number(
                                m_qle_marker_one->text().toULongLong() +
                                m_qle_marker_diff->text().toULongLong()));
}

void zoomToMarker(void)
{
    std::cerr << "WaveWindow::zoomToMarker()" << std::endl;
    m_is_zoom_to_marker = true;
    plot();
}

void findMax(void)
{
    std::vector<double> magnitudes(Dsss::getAbs(m_data));
    uint64_t start = std::min(m_plotter->getMarkerOnePos(),
                              magnitudes.size() - 1);
    uint64_t stop =  std::min(m_plotter->getMarkerTwoPos(),
                              magnitudes.size() - 1);

    if(start > stop)
    {
        uint64_t tmp = start;
        start = stop;
        stop = tmp;
    }

    uint64_t pos =
    std::distance(magnitudes.begin(),
                  std::max_element(magnitudes.begin() + start,
                                   magnitudes.begin() + stop));

    QMessageBox *qmb_info = new QMessageBox;
    qmb_info->setText(QString::number(magnitudes[pos]) + " bei "
                      + QString::number(pos));
    qmb_info->show();
}

public slots:

void
markerWasChanged(void)
{
    m_qle_marker_one->setText(QString::number(m_plotter->getMarkerOnePos()
                                              + m_x_offset));
    m_qle_marker_two->setText(QString::number(m_plotter->getMarkerTwoPos()
                                              + m_x_offset));

    m_qle_marker_diff->setText(QString::number(
                                   m_qle_marker_two->text().toLongLong() -
                                   m_qle_marker_one->text().toLongLong()));
}

/// @brief Setzt internen Datenspeicher
/// @param input Abtastwerte
/// @param leng  Anzahlt der Abtastwerte
void setData(const std::complex<float>* input, uint64_t leng)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::setData(1)" << std::endl;
#endif
    m_data.resize(leng);
    std::copy(input, input + leng, m_data.begin());

    plot();
}
/// @brief Setzt internen Datenspeicher
/// @param input Eingangsvektor mit Abtastwerten
void setData(const std::vector<std::complex<float>> &input)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::setData(2)" << std::endl;
#endif
    m_data = input;
    m_ql_amount_of_samples->setText(QString::number(input.size()));

    plot();
}

/// @brief Setzt die Abtastrate der Darstellung
/// @param samp_rate Abtastrate in Hz
void setSamplerate(double samp_rate)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::setSamplerate()" << std::endl;
#endif
    m_ql_samp_rate->setText(QString::number(samp_rate));

    plot();
}

/// @brief Setzt die Abtastrate der Darstellung
/// @param samp_rate Abtastrate in Hz
void setXOffset(uint64_t offset)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveWindow::setXOffset()" << std::endl;
#endif
    m_x_offset = offset;

    plot();
}

};

#endif // WAVEWINDOW_HPP
