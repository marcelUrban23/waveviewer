#ifndef WAVEPAINTER_HPP
#define WAVEPAINTER_HPP

// Qt
#include <QWidget>
#include <QMouseEvent>
#include <QPainter>
#include <QVector>
#include <QGraphicsView>

// STL
#include <vector>
#include <complex>
#include <iostream>


//#define DEBUG_FUNCTION_CALL

class WavePainter : public QWidget
{
Q_OBJECT

std::vector<double> m_magnitudes;

double  m_x_axes_ratio;

const uint64_t x_height = 50,
               y_width  = 150;

double m_offset,
       m_y_axes_minmax,
       m_min,
       m_max;

uint64_t m_selection_start,
         m_selection_stop;

uint64_t m_x_offset;

int m_marker_one_pos,
    m_marker_two_pos;


QPixmap m_wavemap;




// wurde-veraendert- Indikatoren
bool m_has_marker_changed,
     m_has_wave_changed,

     m_drag_marker_one,
     m_drag_marker_two,

     m_is_selection_area;

int m_marker_status;

QPixmap
getYAxesPixmap(uint64_t width, uint64_t height)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::getYAxesPixmap()" << std::endl;
#endif
    QPixmap axes(width, height);
    axes.fill();
    QPainter painter(&axes);
    QPen pen = painter.pen();
    pen.setWidth(3);
    pen.setColor(Qt::black);
    painter.setPen(pen);
    painter.setOpacity(0.2);
    // x
    painter.drawLine(axes.width(), 0, axes.width(), axes.height());

    return axes;
}


QPixmap
getXAxesPixmap(uint64_t width, uint64_t height)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::getXAxesPixmap()" << std::endl;
#endif
    const int marker = 7;

    QPixmap axes(width, height);
    axes.fill();
    QPainter painter(&axes);
    QPen pen = painter.pen();
    pen.setWidth(3);
    pen.setColor(Qt::black);
    painter.setPen(pen);
    painter.setOpacity(0.2);
    painter.drawLine(0, 0, axes.width(), 0);

    // Werte an die X-Achse anbringen
    for(uint64_t w = 0; w <= marker; ++w)
    {
        // Aequidistante Striche
        painter.drawLine(w * width / marker, 0, w * width / marker, 5);
        // Werte unter den Strichen
        painter.drawText(w * width / marker, 25,
                         QString::number(m_x_offset
                                         + m_magnitudes.size() * w / marker));
    }

    return axes;
}



QPixmap
getWavePixmap(uint64_t width, uint64_t height)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::getWavePixmap()" << std::endl;
#endif
    QPixmap wave(width, height);
    wave.fill();
    QPainter qp_wave(&wave);
    QPen pen = qp_wave.pen();
    pen.setWidth(1);
    pen.setColor(Qt::blue);
    qp_wave.setPen(pen);

    for(uint64_t w = 1; w < m_magnitudes.size(); ++w)
    {
        QPoint x((w - 1) * m_x_axes_ratio,
                 height - (m_magnitudes[w -1] - m_offset) * m_y_axes_minmax);

        QPoint y (w * m_x_axes_ratio,
                  height - (m_magnitudes[w] - m_offset) * m_y_axes_minmax);


        qp_wave.drawLine(x, y);
    }

    return wave;
}

QPixmap
drawMarker(QPixmap wave)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::drawMarker()" << std::endl;
#endif
    QPainter qp(&wave);
    QPen pen = qp.pen();
    pen.setColor(Qt::red);
    qp.setPen(pen);


    qp.drawLine((m_marker_one_pos * m_x_axes_ratio), 0,
                (m_marker_one_pos * m_x_axes_ratio), this->height());
    qp.drawText(m_marker_one_pos * m_x_axes_ratio + 2, 15, "1");

    qp.drawLine((m_marker_two_pos * m_x_axes_ratio), 0,
                (m_marker_two_pos * m_x_axes_ratio), this->height());
    qp.drawText(m_marker_two_pos * m_x_axes_ratio + 2, 15, "2");

    return wave;
}


void
mousePressEvent(QMouseEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::mousePressEvent()" << std::endl;
#endif
    if(event->buttons() & Qt::Key_Shift & Qt::LeftButton)
    {
        m_is_selection_area = true;
    }

    if(event->buttons() & Qt::LeftButton)
    {
        // absolute Amplitude-Position in bestimmen
        int mouse_pos = (event->pos().x() - y_width) / m_x_axes_ratio;

        // Pruefen, ob Marker verschoben werden solln
        if(mouse_pos == m_marker_one_pos)
        {
            m_drag_marker_one = true;
        }
        else if(mouse_pos == m_marker_two_pos)
        {
            m_drag_marker_two = true;
        }
        // setze den Marker
        else
        {
            if(m_marker_status == 0)
            {
                setMarkerOne((event->pos().x() - y_width) / m_x_axes_ratio);
                m_marker_status = 1;
            }
            else if(m_marker_status == 1)
            {
                setMarkerTwo((event->pos().x() - y_width) / m_x_axes_ratio);
                m_marker_status = 0;
            }
        }
    }
}


/// @brief Zustaendig fuer:
///
void
mouseMoveEvent(QMouseEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::mouseMoveEvent()" << std::endl;
#endif

    if(m_magnitudes.empty()) return;

    // Amplituden-Position bestimmen
    int mouse_pos = (event->pos().x() - y_width) / m_x_axes_ratio;
    if(mouse_pos < 0) return;

    if(m_drag_marker_one)
    {
        setMarkerOne((event->pos().x() - y_width) / m_x_axes_ratio);
    }
    if(m_drag_marker_two)
    {
        setMarkerTwo((event->pos().x() - y_width) / m_x_axes_ratio);
    }


    if(mouse_pos == m_marker_one_pos)
    {
        this->setCursor(Qt::SizeHorCursor);
    }
    else if(mouse_pos == m_marker_two_pos)
    {
        this->setCursor(Qt::SizeHorCursor);
    }
    else
    {
        this->setCursor(Qt::ArrowCursor);
    }

    if(event->buttons() & Qt::LeftButton)
    {
        m_selection_start = event->pos().x();
    }

}

void
mouseReleaseEvent(QMouseEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::mouseReleaseEvent()" << std::endl;
#endif
//    if(event->button() & Qt::LeftButton)
//    {
//        m_selection_stop = event->pos().x();
//        QPainter painter(this);
//        painter.drawRect(0, )
//    }
    m_drag_marker_one = false;
    m_drag_marker_two = false;

}



/// @brief Besteht aus drei Bereichen: X- und Y-Achse sowie den Abtastwerten
void
paintEvent(QPaintEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::paintEvent()" << std::endl;
#endif
    if(m_magnitudes.empty()) return;

    QPainter qp_main(this);

    qp_main.fillRect(0, 0, this->width(), this->height(), Qt::white);

    if(m_has_wave_changed)
    {
        m_wavemap = getWavePixmap(this->width() - y_width,
                                  this->height() - x_height);
        m_has_wave_changed = false;
    }

    QPixmap qp_wave_with_marker;
    if(m_has_marker_changed)
    {
        qp_wave_with_marker = drawMarker(m_wavemap);
        emit markerChanged();
        m_has_marker_changed = false;
    }
    else
    {
        qp_wave_with_marker = m_wavemap;
    }

    qp_main.drawPixmap(0, 0, getYAxesPixmap(y_width, this->height() - x_height));
    qp_main.drawPixmap(y_width, this->height() - x_height,
                       getXAxesPixmap(this->width() - y_width, x_height));


    qp_main.drawPixmap(y_width, 0,
                       qp_wave_with_marker);
}


/// @brief Aktualisiert den Veritkalen Offset und die Min-Max Einteilung
void
resizeEvent(QResizeEvent *event)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::resizeEvent()" << std::endl;
#endif

    updateXRatio();
    updateYMinMax();

    m_has_marker_changed = true;
    m_has_wave_changed = true;

    update();
}


void
updateYMinMax(void)
{
    m_y_axes_minmax = (this->height() - x_height) / (m_max - m_min);
    m_offset = m_min;
}

void
updateXRatio(void)
{
    m_x_axes_ratio = static_cast<double>(this->width() - y_width)
                    / static_cast<double>(m_magnitudes.size() - 1);
}

public:

WavePainter(void) :
    m_x_axes_ratio(1.0),
    m_offset(0.0),
    m_y_axes_minmax(1.0),
    m_has_marker_changed(false),
    m_has_wave_changed(false),
    m_drag_marker_one(false),
    m_marker_one_pos(0),
    m_drag_marker_two(false),
    m_marker_two_pos(0),
    m_marker_status(0),
    m_x_offset(0)
{
    QPainter qp_main(this);
    qp_main.fillRect(0, 0, this->width(), this->height(), Qt::white);

    updateXRatio();
    updateYMinMax();

    // Marker-sensitiv
    this->setMouseTracking(true);


}

~WavePainter(void)
{
    m_magnitudes.resize(0);
}


uint64_t getMarkerOnePos(void){return m_marker_one_pos;}
uint64_t getMarkerTwoPos(void){return m_marker_two_pos;}


/// @brief Gibt die dem Pixel naheliegenste Amplitude zrurueck
/// @param pixel relativ zu 0,0
/// @return Amplitude
double
getMagnitudeNearPixel(uint64_t pixel)
{
    uint64_t magnitude_pos = (pixel - y_width) / m_x_axes_ratio;

    if(magnitude_pos < m_magnitudes.size())
    {
        return m_magnitudes[magnitude_pos];
    }

    throw std::runtime_error("FEHLER WavePainter::getMagnitudeNearPixel())");
}

void
setMagnitudes(std::vector<double>::const_iterator input_begin,
              std::vector<double>::const_iterator input_end,
              double x_offset = 0.0, double y_offset = 0.0)
{
#ifdef DEBUG_FUNCTION_CALL
    std::cerr << "WavePainter::setMagnitudes()" << std::endl;
#endif
    m_magnitudes = std::vector<double>(input_begin, input_end);

    m_min = *std::min_element(input_begin, input_end);
    m_max = *std::max_element(input_begin, input_end);

    m_has_wave_changed = true;
//    m_has_marker_changed = true;
    updateYMinMax();
    updateXRatio();

    update();
}


void
setXOffset(uint64_t offset)
{
    m_x_offset = offset;

    update();
}

public slots:


void
addMagnitudes(const std::vector<double> &input)
{
    m_magnitudes.reserve(m_magnitudes.size() + input.size());
    for(auto value : input)
    {
        m_magnitudes.push_back(value);
    }
}

/// @brief Setzt die anzuzeigenden Werte und aktualisiert die x- und y-Ration.
/// @param input Amplituden
void
setMagnitudes(const std::vector<double> &input)
{
    setMagnitudes(input.begin(), input.end());
}

/// @brief
/// @param position
void setMarkerOne(QString position)
{
    setMarkerOne(position.toULongLong());
}
void setMarkerOne(int position)
{
    m_marker_one_pos = std::min(static_cast<int>(m_magnitudes.size()) - 1,
                                position);
    m_has_marker_changed = true;
    update();
}
void setMarkerTwo(QString position)
{
    setMarkerTwo(position.toULongLong());
}
void setMarkerTwo(int position)
{
    m_marker_two_pos = std::min(static_cast<int>(m_magnitudes.size()) - 1,
                                position);
    m_has_marker_changed = true;

    update();
}



signals:

void
sendPointsInArea(const std::vector<double>);

void markerChanged(void);



};


#endif // WAVEPAINTER_HPP
