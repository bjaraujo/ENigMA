
#pragma once

#include <QtWidgets/QGraphicsView>
#include <QtGui/QMouseEvent>
#include <QtGui/QWheelEvent>

class MyGraphicsView : public QGraphicsView
{
    Q_OBJECT
    
public:
    MyGraphicsView();

private:

    double m_s;
    double m_minScale;
    double m_maxScale;

    void scaleBy(qreal scaleFactor);

protected:
    //void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    
};

