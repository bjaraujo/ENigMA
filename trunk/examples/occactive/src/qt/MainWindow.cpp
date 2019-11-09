
#include <QtCore/QElapsedTimer>

#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>
#include <QtGui/QInputDialog>
#include <QtGui/QColorDialog>
#include <QtGui/QMdiSubWindow>

#include "MainWindow.h"

QtMainWindow::QtMainWindow(QApplication* app) : m_shortCutDelete(QShortcut(QKeySequence(Qt::Key_Delete), this, SLOT(editDelete())))
{

    try
    {

        m_app = app;

        this->setupUi(this);

        this->setAttribute(Qt::WA_QuitOnClose);

        this->setWindowTitle(QApplication::applicationName());

        createActions();
        createProgressBar();
        
        readSettings();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

QtMainWindow::~QtMainWindow()
{

    try
    {

        mdiArea->removeSubWindow(m_manager->window());

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::closeEvent(QCloseEvent *event)
{

    try
    {

        writeSettings();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::createActions()
{

    try
    {

        // File
        connect(actionFileNew, SIGNAL(triggered()), this, SLOT(fileNew()));
        connect(actionFileOpen, SIGNAL(triggered()), this, SLOT(fileOpen()));
        connect(actionFileSaveAs, SIGNAL(triggered()), this, SLOT(fileSaveAs()));
        connect(actionFileImport, SIGNAL(triggered()), this, SLOT(fileImport()));
        connect(actionFileExport, SIGNAL(triggered()), this, SLOT(fileExport()));
        connect(actionExit, SIGNAL(triggered()), this, SLOT(applicationExit()));
        
        // Edit
        connect(actionEditDelete, SIGNAL(triggered()), this, SLOT(editDelete()));
        
        // View
        connect(actionViewRenderWireframe, SIGNAL(triggered()), this, SLOT(viewWireframe()));
        connect(actionViewRenderShaded, SIGNAL(triggered()), this, SLOT(viewShaded()));

        connect(actionViewOrientationFront, SIGNAL(triggered()), this, SLOT(viewFront()));
        connect(actionViewOrientationBack, SIGNAL(triggered()), this, SLOT(viewBack()));
        connect(actionViewOrientationLeft, SIGNAL(triggered()), this, SLOT(viewLeft()));
        connect(actionViewOrientationRight, SIGNAL(triggered()), this, SLOT(viewRight()));
        connect(actionViewOrientationTop, SIGNAL(triggered()), this, SLOT(viewTop()));
        connect(actionViewOrientationBottom, SIGNAL(triggered()), this, SLOT(viewBottom()));
        connect(actionViewOrientationIsometric, SIGNAL(triggered()), this, SLOT(viewIsometric()));

        connect(actionViewColor, SIGNAL(triggered()), this, SLOT(viewColor()));

        // Select 
        connect(actionSelectVertices, SIGNAL(triggered()), this, SLOT(selectVertices()));
        connect(actionSelectEdges, SIGNAL(triggered()), this, SLOT(selectEdges()));
        connect(actionSelectWires, SIGNAL(triggered()), this, SLOT(selectWires()));
        connect(actionSelectFaces, SIGNAL(triggered()), this, SLOT(selectFaces()));
        connect(actionSelectSolids, SIGNAL(triggered()), this, SLOT(selectSolids()));
        connect(actionSelectAll, SIGNAL(triggered()), this, SLOT(selectAll()));

        // Shape 
        connect(actionShapeDecompose, SIGNAL(triggered()), this, SLOT(decomposeShape()));
        connect(actionShapeFuse, SIGNAL(triggered()), this, SLOT(fuseShape()));
        connect(actionShapeOffset, SIGNAL(triggered()), this, SLOT(offsetShape()));
        connect(actionShapeBuild, SIGNAL(triggered()), this, SLOT(buildShape()));
        
        // Mesh 
        connect(actionMeshSurface, SIGNAL(triggered()), this, SLOT(meshSurface()));
        connect(actionMeshVolume, SIGNAL(triggered()), this, SLOT(meshVolume()));
        connect(actionMeshStitch, SIGNAL(triggered()), this, SLOT(meshStitch()));

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::createProgressBar()
{

    try
    {

        // Progress bar
        m_progressBar.setVisible(true);
        m_progressBar.setTextVisible(true);

        m_progressBar.setMinimum(0);
        m_progressBar.setMaximum(100);

        infoBar->addPermanentWidget(&m_progressBar);

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::readSettings()
{

    try
    {

        // Read settings INI file
        QSettings settings(QApplication::applicationDirPath().append(QDir::separator()).append("ACTiVE.ini"), QSettings::IniFormat);

        settings.beginGroup("GUI");

        m_curDir = settings.value("CurDir").toString();

        int x = settings.value("PosX").toInt();
        int y = settings.value("PosY").toInt();
        int w = settings.value("Width").toInt();
        int h = settings.value("Height").toInt();

        this->setGeometry(x, y, w, h);
        this->resize(w, h);

        settings.endGroup();

        settings.beginGroup("Mesh");

        m_interval = settings.value("Interval").toInt();
        m_maxElem = settings.value("MaxElem").toInt();

        settings.endGroup();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::writeSettings()
{

    try
    {

        // Write settings INI file
        QSettings settings(QApplication::applicationDirPath().append(QDir::separator()).append("ACTiVE.ini"), QSettings::IniFormat);

        settings.beginGroup("GUI");

        settings.setValue("CurDir", m_curDir);

        settings.setValue("PosX", this->geometry().x());
        settings.setValue("PosY", this->geometry().y());
        settings.setValue("Width", this->geometry().width());
        settings.setValue("Height", this->geometry().height());

        settings.endGroup();

        settings.sync();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::fileNew()
{

    try
    {

        m_manager->reset();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::fileOpen()
{

}

void QtMainWindow::fileSaveAs()
{

}

void QtMainWindow::fileImport()
{

    try
    {
    
        QString selfilter = tr("STEP files (*.stp *.step)");

        QString fileName = QFileDialog::getOpenFileName(this, tr("Import"), m_curDir, 
            tr("All files (*.*);;STEP files (*.stp *.step);;IGES files (*.igs *.iges);;STL files (*.stl)"),
            &selfilter);

        if (fileName == QString::null) return;

        QFileInfo fileInfo(fileName);
        m_curDir = fileInfo.absolutePath();

        if (!QFileInfo(fileName).exists())
        {
            qDebug() << "Error: file does not exist: " << fileName;
            return;
        }

        QApplication::setOverrideCursor(Qt::WaitCursor);            

        QString fileExt = QFileInfo(QString(fileName)).completeSuffix();

        if (QString::compare(fileExt, "stp", Qt::CaseInsensitive) == 0 || QString::compare(fileExt, "step", Qt::CaseInsensitive) == 0)
        {
            // Import STEP file
            m_manager->importGeometry(OccTranslator::FormatSTEP, fileName);
        }
        else if (QString::compare(fileExt, "iges", Qt::CaseInsensitive) == 0 || QString::compare(fileExt, "igs", Qt::CaseInsensitive) == 0)
        {
            // Import IGES file
            m_manager->importGeometry(OccTranslator::FormatIGES, fileName);
        }
        else if (QString::compare(fileExt, "brep", Qt::CaseInsensitive) == 0)
        {
            // Import BREP file
            m_manager->importGeometry(OccTranslator::FormatBREP, fileName);
        }
        else if (QString::compare(fileExt, "csfdb", Qt::CaseInsensitive) == 0)
        {
            // Import CSFDB file
            m_manager->importGeometry(OccTranslator::FormatCSFDB, fileName);
        }
        else if (QString::compare(fileExt, "stl", Qt::CaseInsensitive) == 0)
        {
            // Import STL file
            m_manager->importMesh(OccTranslator::FormatSTL, fileName);
        }
        else
            return;

        QApplication::restoreOverrideCursor();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::fileExport()
{

    try
    {

        QString selfilter = tr("STL files (*.stl)");
        QString fileName = QFileDialog::getSaveFileName(this, tr("Export"), m_curDir, 
            tr("STL files (*.stl);;Gmsh files (*.msh);;Vtk files (*.vtk)"),
            &selfilter);

        if (fileName == QString::null) return;

        QString fileExt = QFileInfo(QString(fileName)).completeSuffix();

        if (QString::compare(fileExt, "stl", Qt::CaseInsensitive) == 0)
        {
            // Export STL
            m_manager->exportMesh(OccTranslator::FormatSTL, fileName);
        }
        else if (QString::compare(fileExt, "msh", Qt::CaseInsensitive) == 0)
        {
            // Export mesh
            m_manager->exportMesh(OccTranslator::FormatGmsh, fileName);
        }
        else if (QString::compare(fileExt, "vtk", Qt::CaseInsensitive) == 0)
        {
            // Export mesh
            m_manager->exportMesh(OccTranslator::FormatVtk, fileName);
        }

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::applicationExit()
{

    try
    {

        m_app->quit();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

// Select

void QtMainWindow::selectVertices()
{

    try
    {

        m_manager->selectVertices();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::selectEdges()
{

    try
    {

        m_manager->selectEdges();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::selectWires()
{

    try
    {

        m_manager->selectWires();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::selectFaces()
{

    try
    {

        m_manager->selectFaces();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::selectSolids()
{

    try
    {

        m_manager->selectSolids();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::selectAll()
{

    try
    {

        m_manager->selectAll();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

// Edit

void QtMainWindow::editDelete()
{

    try
    {

        m_manager->deleteSelected();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

// View

void QtMainWindow::viewWireframe()
{

    try
    {

        m_manager->viewWireframe();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewShaded()
{

    try
    {

        m_manager->viewShaded();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewFront()
{

    try
    {

        m_manager->viewFront();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewBack()
{

    try
    {

        m_manager->viewBack();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewLeft()
{

    try
    {

        m_manager->viewLeft();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewRight()
{

    try
    {

        m_manager->viewRight();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewTop()
{

    try
    {

        m_manager->viewTop();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewBottom()
{

    try
    {

        m_manager->viewBottom();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewIsometric()
{

    try
    {

        m_manager->viewIsometric();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::viewColor()
{

    try
    {

        QColor color = QColorDialog::getColor();
    
        if (color.isValid())
        {    
            m_manager->setColor(color);
        }

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

// Shape
void QtMainWindow::decomposeShape()
{

    try
    {

        m_manager->decomposeShape();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::fuseShape()
{

    try
    {

        m_manager->fuseShape();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::offsetShape()
{

    try
    {

        QString text = QInputDialog::getText(this, "Input", "Enter offset", QLineEdit::Normal, "1.0");

        if (text == QString::null)         
            return;

        m_manager->offsetShape(text.toDouble());

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::buildShape()
{

    try
    {

        m_manager->buildShape();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

// Mesh
void QtMainWindow::meshSurface()
{

    try
    {

        QString text = QInputDialog::getText(this, "Input", "Enter mesh size", QLineEdit::Normal, "2.0");

        if (text == QString::null)
            return;

        m_manager->meshSurface(text.toDouble());

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::meshVolume()
{

    try
    {

        QString text = QInputDialog::getText(this, "Input", "Enter mesh size", QLineEdit::Normal, "2.0");

        if (text == QString::null)
            return;

        m_manager->meshVolume(text.toDouble());

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::meshStitch()
{

    try
    {

        m_manager->meshStitch();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::setManager(QSharedPointer<OccManager> aManager)
{

    try
    {

        m_manager = aManager;

        treeView->setModel(m_manager->treeViewModel());
        treeView->setSelectionModel(m_manager->treeViewSelectionModel());

        connect(m_manager.data(), SIGNAL(onSelectionChanged()), this, SLOT(selectionChanged()));

        // Graphics window
        mdiArea->addSubWindow(m_manager->window());

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

void QtMainWindow::selectionChanged()
{

    try
    {

        this->showProperties();

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << tr("Unknown exception:") << QString(this->metaObject()->className());
    }

}

void QtMainWindow::showProperties()
{

    try
    {

    }
    catch(const std::exception & ex)
    {
        qDebug() << QString(ex.what());
    }
    catch(...)
    {
        qDebug() << "Unknown exception: " + QString(this->metaObject()->className());
    }

}

