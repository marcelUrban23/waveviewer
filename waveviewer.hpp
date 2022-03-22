#ifndef WAVEVIEWER_HPP
#define WAVEVIEWER_HPP

#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMimeData>

#include "wavedata.hpp"
#include "wavewindow.hpp"
#include "filterandresampler.hpp"


//#define DEBUG_FUNCTION_CALL


/// @author Marcel Urban
/// @brief Diese Klasse stellt den aeusseren Rahmen dar.
///        Hier werden folgende Dinge gehandhabt:
///        Quelle
///        Anzahl und Orientierung der Abtastwerte
///        Datentyp
///        Byte-Offset der Quelle (Datei)
///        Font-Groesse
class WaveViewer : public QWidget
{
Q_OBJECT

QVBoxLayout *m_qvbl_main,
            *m_qvbl_viewer;

WaveData    *m_wavedata;


FilterAndResampler *m_filter_and_resampler;

QDataStream *m_in;

QFile *m_file;

QLineEdit *m_qle_sample_rate,
          *m_qle_sample_all,
          *m_qle_sample_amount,
          *m_qle_sample_begin,
          *m_qle_sample_end,
          *m_qle_filepath,
          *m_qle_byteoffset,

// Filtern und Resampeln
          *m_qle_shift,
          *m_qle_transition,
          *m_qle_passband,
          *m_qle_resample_frq;
QCheckBox *m_qcb_enable_shift,
          *m_qcb_enable_filter,
          *m_qcb_enable_resample;

QString m_file_name;

QComboBox *m_qcb_datatype;

QScrollBar *m_hz_scrollbar;

bool m_is_selection_area,
     m_is_file;

uint64_t m_fontsize;

public:

WaveViewer(void) :
    m_is_file(false),
    m_fontsize(7)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::WaveViewer()" << std::endl;
#endif
    setGeometry(0, 0, 500, 900);

    m_filter_and_resampler = new FilterAndResampler();

    m_wavedata = new WaveData;

    createGUI();

    setAcceptDrops(true);

    setNewFont();
}
~WaveViewer()
{

}

private:

void
createGUI(void)
{
    QPushButton *qpb_open_file = new QPushButton;
    qpb_open_file->setText("Datei");
    connect(qpb_open_file, SIGNAL(clicked()),
            this         , SLOT(setFile()));

    QPushButton *qpb_new_wavewindow = new QPushButton("neues Fenster");
    connect(qpb_new_wavewindow, SIGNAL(clicked()),
            this              , SLOT(newWavewindow()));

    m_qle_filepath = new QLineEdit;
    connect(m_qle_filepath, SIGNAL(editingFinished()),
            this          , SLOT(setFilename()));

    m_qle_sample_rate = new QLineEdit;
    m_qle_sample_rate->setFixedWidth(100);
    m_qle_sample_rate->setText("1000000");
    connect(m_qle_sample_rate, SIGNAL(editingFinished()),
            this             , SLOT(setSamplerate()));

    m_qcb_datatype = new QComboBox;
    for(auto data_type : m_wavedata->getDatatypes())
    {
        m_qcb_datatype->addItem(QString::fromStdString(data_type));
    }
    connect(m_qcb_datatype, SIGNAL(currentIndexChanged(int)),
            this           , SLOT(setDataType(int)));


    m_qle_byteoffset = new QLineEdit;
    m_qle_byteoffset->setFixedWidth(30);
    m_qle_byteoffset->setText("0");
    connect(m_qle_byteoffset, SIGNAL(editingFinished()),
            this            , SLOT(setByteoffset()));

    QCheckBox *qcb_swap_bytes = new QCheckBox;
    connect(qcb_swap_bytes, SIGNAL(clicked(bool)),
            this          , SLOT(setSwappBytes(bool)));

    QHBoxLayout *qhbl_file = new QHBoxLayout;
    qhbl_file->addWidget(qpb_open_file);
    qhbl_file->addWidget(m_qle_filepath);
    qhbl_file->addWidget(new QLabel("Abtastrate: "));
    qhbl_file->addWidget(m_qle_sample_rate);
    qhbl_file->addWidget(new QLabel("Datentyp: "));
    qhbl_file->addWidget(m_qcb_datatype);
    qhbl_file->addWidget(new QLabel("Byteoffset: "));
    qhbl_file->addWidget(m_qle_byteoffset);
    qhbl_file->addWidget(new QLabel("Swap Bytes: "));
    qhbl_file->addWidget(qcb_swap_bytes);


    m_qle_sample_all = new QLineEdit;
    m_qle_sample_all->setFixedWidth(100);
    m_qle_sample_all->setReadOnly(true);

    m_qle_sample_begin = new QLineEdit;
    m_qle_sample_begin->setFixedWidth(100);
    m_qle_sample_begin->setText("0");
    connect(m_qle_sample_begin, SIGNAL(editingFinished()),
            this             , SLOT(setSamplerange()));

    m_qle_sample_end   = new QLineEdit;
    m_qle_sample_end->setFixedWidth(100);
    m_qle_sample_end->setReadOnly(true);
    m_qle_sample_end->setText("8191");

    m_qle_sample_amount = new QLineEdit;
    m_qle_sample_amount->setFixedWidth(100);
    m_qle_sample_amount->setText("8192");
    connect(m_qle_sample_amount, SIGNAL(editingFinished()),
            this               , SLOT(setSamplerange()));

    QHBoxLayout *qhbl_data_read_settings = new QHBoxLayout;
    qhbl_data_read_settings->addWidget(qpb_new_wavewindow);
    qhbl_data_read_settings->addWidget(new QLabel("Von: "));
    qhbl_data_read_settings->addWidget(m_qle_sample_begin);
    qhbl_data_read_settings->addWidget(new QLabel("Anzahl: "));
    qhbl_data_read_settings->addWidget(m_qle_sample_amount);
    qhbl_data_read_settings->addWidget(new QLabel("Bis: "));
    qhbl_data_read_settings->addWidget(m_qle_sample_end);
    qhbl_data_read_settings->addWidget(new QLabel("Gesamt: "));
    qhbl_data_read_settings->addWidget(m_qle_sample_all);
    qhbl_data_read_settings->addStretch();

    
    
    QPushButton *qpb_save_data_to_file = new QPushButton("speichern als...");
//    connect(qpb_save_data_to_file, SIGNAL(clicked(bool)),
//            this       , SLOT(zoom()));
    
    QHBoxLayout *qhbl_data_write_to_file = new QHBoxLayout;
    qhbl_data_write_to_file->addWidget(qpb_save_data_to_file);
    
    m_qcb_enable_shift = new QCheckBox;
    m_qcb_enable_shift->setChecked(false);
    connect(m_qcb_enable_shift, SIGNAL(stateChanged(int)),
            this              , SLOT(sendData()));

    m_qle_shift = new QLineEdit;
    m_qle_shift->setFixedWidth(100);
    m_qle_shift->setText("0.0");
    connect(m_qle_shift, SIGNAL(editingFinished()),
            this              , SLOT(sendData()));

    m_qcb_enable_filter = new QCheckBox;
    m_qcb_enable_filter->setChecked(false);
    connect(m_qcb_enable_filter, SIGNAL(stateChanged(int)),
            this              , SLOT(sendData()));

    m_qle_transition = new QLineEdit;
    m_qle_transition->setFixedWidth(100);
    m_qle_transition->setText("0.0");
    connect(m_qle_shift, SIGNAL(editingFinished()),
            this              , SLOT(sendData()));
    m_qle_passband = new QLineEdit;
    m_qle_passband->setFixedWidth(100);
    m_qle_passband->setText("0.0");
    connect(m_qle_passband, SIGNAL(editingFinished()),
            this              , SLOT(sendData()));

    m_qcb_enable_resample = new QCheckBox;
    m_qcb_enable_resample->setChecked(false);
    connect(m_qcb_enable_resample, SIGNAL(stateChanged(int)),
            this              , SLOT(sendData()));

    m_qle_resample_frq = new QLineEdit;
    m_qle_resample_frq->setFixedWidth(100);
    m_qle_resample_frq->setText("0.0");
    connect(m_qle_resample_frq, SIGNAL(editingFinished()),
            this              , SLOT(sendData()));

    QHBoxLayout *qhbl_data_modify_settings = new QHBoxLayout;
    qhbl_data_modify_settings->addWidget(m_qcb_enable_shift);
    qhbl_data_modify_settings->addWidget(new QLabel("shift: "));
    qhbl_data_modify_settings->addWidget(m_qle_shift);
    qhbl_data_modify_settings->addWidget(m_qcb_enable_filter);
    qhbl_data_modify_settings->addWidget(new QLabel("transition: "));
    qhbl_data_modify_settings->addWidget(m_qle_transition);
    qhbl_data_modify_settings->addWidget(new QLabel("passband: "));
    qhbl_data_modify_settings->addWidget(m_qle_passband);
    qhbl_data_modify_settings->addWidget(m_qcb_enable_resample);
    qhbl_data_modify_settings->addWidget(new QLabel("resample Frq: "));
    qhbl_data_modify_settings->addWidget(m_qle_resample_frq);
    qhbl_data_modify_settings->addStretch();
    

    QVBoxLayout *qvbl_settings = new QVBoxLayout;
    qvbl_settings->addLayout(qhbl_data_read_settings);
    qvbl_settings->addLayout(qhbl_data_modify_settings);

    m_qvbl_viewer = new QVBoxLayout;
    newWavewindow();

    m_hz_scrollbar = new QScrollBar(Qt::Orientation::Horizontal);
    m_hz_scrollbar->setRange(0, 0);
    connect(m_hz_scrollbar, SIGNAL(valueChanged(int)),
            this          , SLOT(handleScrollbar(int)));

    m_qvbl_main = new QVBoxLayout;
    m_qvbl_main->addLayout(qhbl_file);
    m_qvbl_main->addLayout(qvbl_settings);
    m_qvbl_main->addLayout(m_qvbl_viewer);
    m_qvbl_main->addWidget(m_hz_scrollbar);

    setLayout(m_qvbl_main);
}

void
keyPressEvent(QKeyEvent * event)
{
    if(event->key() == 43)
    {
        ++m_fontsize;
        setNewFont();
    }
    else if(event->key() == 45)
    {
        --m_fontsize;
        setNewFont();
    }
}

void
setNewFont(void)
{
    if(m_fontsize > 40) m_fontsize = 40;
    if(m_fontsize <  5) m_fontsize =  5;

    QFont qf_font;
    qf_font.setPointSize(m_fontsize);
    this->setFont(qf_font);
}

private slots:


/// @brief Reagiert auf Veraenderung der Scrollbar
/// @param pos aktuelle Position der Scrollbar
void
handleScrollbar(int pos)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::handleScrollbar()" << std::endl;
#endif
    m_qle_sample_begin->setText(QString::number(pos));
    m_qle_sample_end->setText(QString::number(pos +
                                              m_qle_sample_amount->text().toInt()
                                              - 1));
    sendData();
}


/// @brief Reagiert auf Aenderung des Datentypes und fordert entsprechend neue
///        Daten an.
void
setDataType(int index)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "MainWindow::setDataType()" << std::endl;
#endif
    if( ! m_is_file){return;}
    m_wavedata->setDataType(static_cast<uint32_t>(index));
    m_qcb_datatype->setCurrentIndex(static_cast<int>(index));
    emit pushData(m_wavedata->getSamples(m_qle_sample_begin->text().toInt(),
                                         m_qle_sample_amount->text().toInt()));

}


/// @brief Reagiert auf Reihenfolge der Bytes und fordert entsprechend neue
///        Daten an.
void
setSwappBytes(bool state)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "MainWindow::setSwappBytes()" << std::endl;
#endif
    if( ! m_is_file){return;}
    m_wavedata->setSwapped(state);
    sendData();
}


/// @brief Reagiert auf Anhebung des Startbytes und fordert entsprechend neue
///        Daten an.
void
setByteoffset(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "MainWindow::setByteoffset()" << std::endl;
#endif
    if( ! m_is_file){return;}
    m_wavedata->setByteOffset(m_qle_byteoffset->text().toInt());
    sendData();
}


void
sendData(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::sendData()" << std::endl;
#endif
    std::vector<std::complex<float>> data(
                m_wavedata->getSamples(m_qle_sample_begin->text().toInt(),
                                       m_qle_sample_amount->text().toInt()));

    if(m_qcb_enable_shift->isChecked())
    {
        data = m_filter_and_resampler->shiftFrequency(data,
                                        m_qle_sample_rate->text().toDouble(),
                                        m_qle_shift->text().toDouble());
    }

    if(m_qcb_enable_filter->isChecked())
    {
        data = m_filter_and_resampler->getFilteredSignal(data,
                                           m_qle_passband->text().toDouble(),
                                           m_qle_transition->text().toDouble());
    }

    if(m_qcb_enable_resample->isChecked())
    {
        data = m_filter_and_resampler->resampleCubic(data,
                                         m_qle_sample_rate->text().toDouble()
                                       / m_qle_resample_frq->text().toDouble());
    }
    emit pushXOffset(m_qle_sample_begin->text().toULongLong());
    emit pushData(data);
}

/// @brief Reagiert auf die Anforderung eines neuen Fensters.
void
newWavewindow(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::newWavewindow()" << std::endl;
#endif
    WaveWindow *wv = new WaveWindow;
    m_qvbl_viewer->addWidget(wv);
    connect(this, SIGNAL(pushData(std::vector<std::complex<float>>)),
            wv  , SLOT(setData(std::vector<std::complex<float>>)));

    connect(this, SIGNAL(pushXOffset(uint64_t)),
            wv  , SLOT(setXOffset(uint64_t)));

    connect(this, SIGNAL(pushSamplerate(double)),
            wv  , SLOT(setSamplerate(double)));

    if( ! m_is_file){return;}
    sendData();
    setSamplerate();
}


/// @brief Setzt die Laenge der Scrollbar relativ zur angezeigten bzw. absoluten
///        Anzahl an Abtastwerten der Datei.
void
setScrollbar(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::setScrollbar()" << std::endl;
#endif
    m_hz_scrollbar->setRange(0, m_wavedata->getFileSizeInSamples());
}


/// @brief Reagiert auf Aenderung des Dateipfades in der Adressleiste und
///        versucht neue Datei zu oeffnen bzw. einzulesen.
void
setFilename(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::setFile()" << std::endl;
#endif
    setFile(m_qle_filepath->text());
}


void
setFile(QString filename)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::setFile(1)" << std::endl;
#endif

    if(m_is_file){delete m_wavedata;}

    m_file_name = filename;
    m_wavedata = new WaveData(m_file_name.toStdString());
    if( ! m_wavedata->getIsInputFile())
    {
        QMessageBox info;
        info.setText("keine Datei gefunden");
        info.show();
        return;
    }
    m_is_file = true;

    m_qle_filepath->setText(filename);
    setDataType(m_wavedata->getDataTypeIndex());
    m_qle_sample_all->setText(QString::number(m_wavedata->getFileSizeInSamples()));

    m_qle_sample_rate->setText(QString::number(m_wavedata->m_samp_rate));
    setSamplerate();
    sendData();
    setScrollbar();
}
void
setFile(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::setFile(2)" << std::endl;
#endif

    QString filename = QFileDialog::getOpenFileName(this, ("Open File"),
                                                    "/home/Schreibtisch",
                                                    ("*"));
    setFile(filename);

    this->setWindowTitle(filename);
}



/// @brief Reagiert auf die Aenderung der anzugeigenden Anzahl an Abtastwerten
void
setSamplerange(void)
{
    if( ! m_is_file){return;}
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::setSamplerange()" << std::endl;
#endif

    m_qle_sample_end->setText(QString::number(m_qle_sample_begin->text().toInt()+
                                              m_qle_sample_amount->text().toInt()
                                              - 1));
    sendData();
}


/// @brief Reagiert auf manuelle Aenderung der Abtastrate, dies beinflusst
///        lediglich die Frequenzskala
void
setSamplerate(void)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::setSamplerate()" << std::endl;
#endif

    emit pushSamplerate(m_qle_sample_rate->text().toDouble());
}


/// @brief Reagiert, sobald eine Datei innerhalb des Fensters gezogen wird
void
dragEnterEvent(QDragEnterEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::dragEnterEvent()" << std::endl;
#endif
    event->acceptProposedAction();
}


/// @brief Reagiert, sobald eine Datei innerhalb des Fensters losgelassen wird.
///        Es wird versucht, diese Datei zu oeffnen bzw. zu inzulesen.
void
dropEvent(QDropEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WaveViewer::dropEvent()" << std::endl;
#endif
    foreach (const QUrl &url, event->mimeData()->urls())
    {
        setFile(url.toLocalFile());
    }

}


signals:
    void pushData(std::vector<std::complex<float>> input);

    void pushXOffset(uint64_t offset);

    void pushSamplerate(double samp_rate);

};





#endif // WAVEVIEWER_HPP
