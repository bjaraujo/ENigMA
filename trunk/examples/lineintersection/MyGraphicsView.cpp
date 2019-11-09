// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <cmath>

#include "MyGraphicsView.h"

MyGraphicsView::MyGraphicsView()
{

    m_s = 1.0;
 
    m_minScale = 0.01;
    m_maxScale = 100.0;

    //this->setRenderHints(QPainter::Antialiasing);
    //this->centerOn(0.0, 0.0);

    this->setDragMode(ScrollHandDrag);

}

/*
void MyGraphicsView::mouseMoveEvent(QMouseEvent *event)
{

    QGraphicsView::mouseMoveEvent(event);
 
}
*/

void MyGraphicsView::wheelEvent(QWheelEvent *event)
{

    QGraphicsView::wheelEvent(event);

    scaleBy(std::pow(4.0 / 3.0, (-event->delta() / 240.0)));

}

void MyGraphicsView::scaleBy(qreal scaleFactor)
{

    qreal curScaleFactor = transform().m11();
    if (((curScaleFactor == m_minScale) && (scaleFactor < 1.0)) ||
        ((curScaleFactor == m_maxScale) && (scaleFactor > 1.0))) 
        return;
 
    qreal sc = scaleFactor;
 
    if ((curScaleFactor * sc < m_minScale)&&(sc < 1.0))
    {
        sc = m_minScale / curScaleFactor;
    }
    else
        if ((curScaleFactor * sc > m_maxScale)&&(sc > 1.0))
        {
        sc = m_maxScale / curScaleFactor;
    }

    this->scale(sc, sc);

}
