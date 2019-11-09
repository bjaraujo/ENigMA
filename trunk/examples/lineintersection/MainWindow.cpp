// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>

#include <QtGui>

#include "MainWindow.h"
#include "MyGraphicsView.h"

#include "GeoLine.hpp"

using namespace ENigMA::geometry;

MainWindow::MainWindow(QWidget *parent)
{

    Ui::MainWindow::setupUi(this);

    MyGraphicsView* graphicsView = new MyGraphicsView;

    this->setCentralWidget(graphicsView);

    QGraphicsScene* scene = new QGraphicsScene;
    scene->setSceneRect(-50.0, -50.0, 50.0, 50.0);

    graphicsView->setScene(scene);

    CGeoCoordinate<float> aVertex1(2.69535, 1.0938, 0);
    CGeoCoordinate<float> aVertex2(3.5, 0.606218, 0);

    CGeoCoordinate<float> aVertex3(3.2, 1.03781, 0);
    CGeoCoordinate<float> aVertex4(2.5, 0.606218, 0);

    // Qt interection

    QLineF aQtLine1;
    aQtLine1.setP1(QPointF(aVertex1.x(), aVertex1.y()));
    aQtLine1.setP2(QPointF(aVertex2.x(), aVertex2.y()));

    QLineF aQtLine2;
    aQtLine2.setP1(QPointF(aVertex3.x(), aVertex3.y()));
    aQtLine2.setP2(QPointF(aVertex4.x(), aVertex4.y()));

    scene->addLine(aQtLine1, QPen(Qt::red));
    scene->addLine(aQtLine2, QPen(Qt::green));

    QPointF aQtPoint;

    if (aQtLine1.intersect(aQtLine2, &aQtPoint) == QLineF::BoundedIntersection)
    {

        std::cout << "Qt: " << aQtPoint.x() << "," << aQtPoint.y() << std::endl;

        double rad = 1;
        scene->addEllipse(aQtPoint.x() - rad, aQtPoint.y() - rad, rad*2.0, rad*2.0, QPen(Qt::blue), QBrush(Qt::blue, Qt::SolidPattern));
    }
    else
    {
        std::cout << "Do not intersect!" << std::endl;
    }

    // My intersection    
    CGeoLine<float> aLine1;
    aLine1.setStartPoint(aVertex1);
    aLine1.setEndPoint(aVertex2);

    CGeoLine<float> aLine2;
    aLine2.setStartPoint(aVertex3);
    aLine2.setEndPoint(aVertex4);

    CGeoCoordinate<float> aPoint;

    CGeoIntersectionType anIntersectionType;

    if (aLine1.intersects(aLine2, aPoint, anIntersectionType, 1E-4))
    {

        std::cout << "ENigMA: " << aPoint.x() << "," << aPoint.y() << std::endl;

        double rad = 1;
        scene->addEllipse(aPoint.x() - rad, aPoint.y() - rad, rad*2.0, rad*2.0, QPen(Qt::black), QBrush(Qt::black, Qt::SolidPattern));

        if (anIntersectionType == IT_VERTEX)
            std::cout << "--> IT_VERTEX" << std::endl;
        else if (anIntersectionType == IT_EDGE)
            std::cout << "--> IT_EDGE" << std::endl;
        else if (anIntersectionType == IT_COINCIDENT)
            std::cout << "--> IT_COINCIDENT" << std::endl;
        else if (anIntersectionType == IT_INTERNAL)
            std::cout << "--> IT_INTERNAL" << std::endl;
        else if (anIntersectionType == IT_SWAP)
            std::cout << "--> IT_SWAP" << std::endl;
        else if (anIntersectionType == IT_NONE)
            std::cout << "--> IT_NONE" << std::endl;

    }
    else
    {
        std::cout << "Do not intersect!" << std::endl;
    }

}

