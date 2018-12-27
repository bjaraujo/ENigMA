#if !defined WNT
#define QT_CLEAN_NAMESPACE         /* avoid definition of INT32 and INT8 */
#endif

#include <QApplication>
#include <QPainter>
#include <QMenu>
#include <QColorDialog>
#include <QCursor>
#include <QFileInfo>
#include <QMouseEvent>
#include <QRubberBand>

#include <Visual3d_View.hxx>
#include <Graphic3d_ExportFormat.hxx>
#include <Graphic3d_GraphicDriver.hxx>

#include <QWindowsStyle>

#if defined(_WIN32) || defined(__WIN32__)
#include <WNT_Window.hxx>
#elif defined(__APPLE__) && !defined(MACOSX_USE_GLX)
#include <Cocoa_Window.hxx>
#else
#include <QX11Info>
#include <GL/glx.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/Xmu/StdCmap.h>
#include <X11/Xlib.h>
#include <Xw_Window.hxx>
#include <QColormap>
#endif

#include <Aspect_DisplayConnection.hxx>

#include "OccWidget.h"

// the key for multi selection :
#define MULTISELECTIONKEY Qt::ShiftModifier

// the key for shortcut (use to activate dynamic rotation, panning)
#define CASCADESHORTCUTKEY Qt::ControlModifier

// for elastic bean selection
#define ValZWMin 1

static QCursor* defCursor     = NULL;
static QCursor* handCursor    = NULL;
static QCursor* panCursor     = NULL;
static QCursor* globPanCursor = NULL;
static QCursor* zoomCursor    = NULL;
static QCursor* rotCursor     = NULL;

OccWidget::OccWidget() : m_viewActions(0)
{



}

OccWidget::~OccWidget()
{


}

void OccWidget::setContext(Handle(AIS_InteractiveContext) theContext)
{

#if !defined(_WIN32) && !defined(__WIN32__) && (!defined(__APPLE__) || defined(MACOSX_USE_GLX))
    //XSynchronize(x11Display(),true); // it is possible to use QApplication::syncX();
    XSynchronize(x11Info().display(),true); // it is possible to use QApplication::syncX();
#endif
    m_first = true;
    m_context = theContext;

    m_xmin = 0;
    m_ymin = 0;
    m_xmax = 0;
    m_ymax = 0;
    m_curZoom = 0;
    m_rectBand = 0;

    setAttribute(Qt::WA_PaintOnScreen);
    setAttribute(Qt::WA_NoSystemBackground);

#if !defined(_WIN32) && !defined(__WIN32__) && (!defined(__APPLE__) || defined(MACOSX_USE_GLX))
    XVisualInfo* pVisualInfo;
    if (x11Info().display())
    {
        /* Initialization with the default VisualID */
        Visual *v = DefaultVisual(x11Info().display(), DefaultScreen(x11Info().display()));
        int visualID = XVisualIDFromVisual(v);

        /*  Here we use the settings from Optimizer_ViewInfo::TxglCreateWindow() */
        int visualAttr[] = { GLX_RGBA, GLX_DEPTH_SIZE, 1, GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1,
            GLX_BLUE_SIZE, 1, GLX_DOUBLEBUFFER, None };
        pVisualInfo = ::glXChooseVisual(x11Info().display(), DefaultScreen(x11Info().display()), visualAttr);

        if (isVisible())
            hide();

        XSetWindowAttributes a;

        Window p = RootWindow(x11Info().display(), DefaultScreen(x11Info().display()));
        a.colormap = XCreateColormap(x11Info().display(), RootWindow(x11Info().display(), pVisualInfo->screen),
            pVisualInfo->visual, AllocNone);

        QColor color = palette().color(backgroundRole());
        QColormap colmap = QColormap::instance();
        a.background_pixel = colmap.pixel(color);
        a.border_pixel = colmap.pixel(Qt::black);
        if (parentWidget())
            p = parentWidget()->winId();

        Window w = XCreateWindow(x11Info().display(), p,  x(), y(), width(), height(),
            0, pVisualInfo->depth, InputOutput,  pVisualInfo->visual,
            CWBackPixel | CWBorderPixel | CWColormap, &a);

        Window *cmw;
        Window *cmwret;
        int count;
        if (XGetWMColormapWindows(x11Info().display(), topLevelWidget()->winId(), &cmwret, &count))
        {
            cmw = new Window[count+1];
            memcpy((char *)cmw, (char *)cmwret, sizeof(Window)*count);
            XFree((char *)cmwret);
            int i;
            for (i = 0; i < count; ++i)
            {
                if (cmw[i] == winId())  /* replace old window */
                {
                    cmw[i] = w;
                    break;
                }
            }
            if (i >= count)             /* append new window */
                cmw[count++] = w;
        }
        else
        {
            count = 1;
            cmw = new Window[count];
            cmw[0] = w;
        }
        /* Creating new window (with good VisualID) for this widget */
        create(w);
        XSetWMColormapWindows(x11Info().display(), topLevelWidget()->winId(), cmw, count);
        delete [] cmw;

        if (isVisible())
            show();
        if (pVisualInfo)
            XFree((char *)pVisualInfo);
        XFlush(x11Info().display());
    }
#endif

    m_currentMode = CurAction3d_Nothing;
    m_hlrModeIsOn = Standard_False;
    setMouseTracking(true);

    initViewActions();
    initCursors();

    setBackgroundRole(QPalette::NoRole);//NoBackground);
    // set focus policy to threat QContextMenuEvent from keyboard  
    setFocusPolicy(Qt::StrongFocus);
    setAttribute(Qt::WA_PaintOnScreen);
    setAttribute(Qt::WA_NoSystemBackground);

}

void OccWidget::init()
{

    if (m_view.IsNull())
        m_view = m_context->CurrentViewer()->CreateView();

#if defined(_WIN32) || defined(__WIN32__)
    Aspect_Handle aWindowHandle = (Aspect_Handle)winId();
    Handle(WNT_Window) hWnd = new WNT_Window (aWindowHandle);
#elif defined(__APPLE__) && !defined(MACOSX_USE_GLX)
    NSView* aViewHandle = (NSView*)winId();
    Handle(Cocoa_Window) hWnd = new Cocoa_Window (aViewHandle);
#else
    Window aWindowHandle = (Window)winId();
    Handle(Aspect_DisplayConnection) aDispConnection = m_context->CurrentViewer()->Driver()->GetDisplayConnection();
    Handle(Xw_Window) hWnd = new Xw_Window (aDispConnection, aWindowHandle);
#endif // WNT
    m_view->SetWindow(hWnd);
    if (!hWnd->IsMapped())
    {
        hWnd->Map();
    }

    m_view->SetBgGradientColors(Quantity_Color(0.1, 0.1, 0.2, Quantity_TOC_RGB), Quantity_Color(0.7, 0.7, 0.8, Quantity_TOC_RGB), Aspect_GFM_VER, Standard_False);
    m_view->MustBeResized();

    // Axis in the corner
    m_view->ZBufferTriedronSetup();
    m_view->TriedronDisplay(Aspect_TOTP_RIGHT_LOWER, Quantity_NOC_BLACK, 0.1, V3d_ZBUFFER);

}

void OccWidget::paintEvent(QPaintEvent * )
{

    //  QApplication::syncX();
    if(m_first)
    {
        init();
        m_first = false;
    }
    m_view->Redraw();

}

void OccWidget::resizeEvent(QResizeEvent *)
{
    //  QApplication::syncX();
    if(!m_view.IsNull())
    {
        m_view->MustBeResized();
    }
}

void OccWidget::fitAll()
{
    m_view->FitAll();
    m_view->ZFitAll();
    m_view->Redraw();
}

void OccWidget::fitArea()
{

    m_currentMode = CurAction3d_WindowZooming;
    activateCursor(m_currentMode);

}

void OccWidget::zoom()
{

    m_currentMode = CurAction3d_DynamicZooming;
    activateCursor(m_currentMode);

}

void OccWidget::pan()
{

    m_currentMode = CurAction3d_DynamicPanning;
    activateCursor(m_currentMode);

}

void OccWidget::rotation()
{
    
    m_currentMode = CurAction3d_DynamicRotation;
    activateCursor(m_currentMode);

}

void OccWidget::globalPan()
{

    // save the current zoom value
    m_curZoom = m_view->Scale();
    // Do a Global Zoom
    m_view->FitAll();
    // Set the mode
    m_currentMode = CurAction3d_GlobalPanning;
    activateCursor(m_currentMode);

}

void OccWidget::front()
{

    m_view->SetProj(V3d_Xpos);

}

void OccWidget::back()
{

    m_view->SetProj(V3d_Xneg);

}

void OccWidget::top()
{

    m_view->SetProj(V3d_Zpos);

}

void OccWidget::bottom()
{

    m_view->SetProj(V3d_Zneg);

}

void OccWidget::left()
{

    m_view->SetProj(V3d_Ypos);

}

void OccWidget::right()
{

    m_view->SetProj(V3d_Yneg);

}

void OccWidget::axo()
{

    m_view->SetProj(V3d_XposYnegZpos);

}

void OccWidget::reset()
{

    m_view->Reset();

}

void OccWidget::hlrOff()
{

    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_hlrModeIsOn = Standard_False;
    m_view->SetComputedMode (m_hlrModeIsOn);
    QApplication::restoreOverrideCursor();

}

void OccWidget::hlrOn()
{

    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_hlrModeIsOn = Standard_True;
    m_view->SetComputedMode (m_hlrModeIsOn);
    QApplication::restoreOverrideCursor();

}

void OccWidget::updateToggled(bool isOn)
{

    QAction* sentBy = (QAction*)sender();

    if(!isOn)
        return;

    for (int i = ViewFitAllId; i < ViewHlrOffId; ++i)
    {
        QAction* anAction = m_viewActions->at(i);
        if ((anAction == m_viewActions->at(ViewFitAreaId)) ||
            (anAction == m_viewActions->at(ViewZoomId)) ||
            (anAction == m_viewActions->at(ViewPanId)) ||
            (anAction == m_viewActions->at(ViewGlobalPanId)) ||
            (anAction == m_viewActions->at(ViewRotationId)))
        {
            if (anAction && (anAction != sentBy))
            {
                anAction->setCheckable(true);
                anAction->setChecked(false);
            }
            else
            {
                if (sentBy == m_viewActions->at(ViewFitAreaId))
                    setCursor(*handCursor);
                else if    (sentBy == m_viewActions->at(ViewZoomId))
                    setCursor(*zoomCursor);
                else if    (sentBy == m_viewActions->at(ViewPanId))
                    setCursor(*panCursor);
                else if    (sentBy == m_viewActions->at(ViewGlobalPanId))
                    setCursor(*globPanCursor);
                else if (sentBy == m_viewActions->at(ViewRotationId))
                    setCursor(*rotCursor);
                else
                    setCursor(*defCursor);

                sentBy->setCheckable(false);
            }
        }
    }
}

void OccWidget::initCursors()
{

    if (!defCursor)
        defCursor = new QCursor(Qt::ArrowCursor);

    if (!handCursor)
        handCursor = new QCursor(Qt::PointingHandCursor);

    if (!panCursor)
        panCursor = new QCursor(Qt::SizeAllCursor);

    if (!globPanCursor)
        globPanCursor = new QCursor(Qt::CrossCursor);
    
    /*
    if (!zoomCursor)
        zoomCursor = new QCursor(QPixmap(ApplicationCommonWindow::getResourceDir() + QString("/") + QObject::tr("ICON_CURSOR_ZOOM")));
    
    if (!rotCursor)
        rotCursor = new QCursor(QPixmap(ApplicationCommonWindow::getResourceDir() + QString("/") + QObject::tr("ICON_CURSOR_ROTATE")));
    */

}

QList<QAction*>* OccWidget::getViewActions()
{

    initViewActions();
    return m_viewActions;

}

/*!
Get paint engine for the OpenGL viewer. [ virtual public ]
*/
QPaintEngine* OccWidget::paintEngine() const
{

    return 0;

}

void OccWidget::initViewActions()
{

    if (m_viewActions)
        return;

}

void OccWidget::mousePressEvent(QMouseEvent* e)
{

    if (e->button() == Qt::LeftButton)
        onLButtonDown((e->buttons() | e->modifiers()), e->pos());
    else if (e->button() == Qt::MidButton)
        onMButtonDown(e->buttons() | e->modifiers(), e->pos());
    else if (e->button() == Qt::RightButton)
        onRButtonDown(e->buttons() | e->modifiers(), e->pos());

}

void OccWidget::mouseReleaseEvent(QMouseEvent* e)
{

    if (e->button() == Qt::LeftButton)
        onLButtonUp(e->buttons(), e->pos());
    else if (e->button() == Qt::MidButton)
        onMButtonUp(e->buttons(), e->pos());
    else if(e->button() == Qt::RightButton)
        onRButtonUp(e->buttons(), e->pos());

}

void OccWidget::mouseMoveEvent(QMouseEvent* e)
{

    onMouseMove(e->buttons(), e->pos());

}

void OccWidget::wheelEvent(QWheelEvent *e)
{
    
    m_view->Zoom(0.0, 0.0, e->delta() * 0.1, 0.0);

}

void OccWidget::activateCursor(const CurrentAction3d mode)
{

    switch(mode)
    {
    case CurAction3d_DynamicPanning:
        setCursor(*panCursor);
        break;
    case CurAction3d_DynamicZooming:
        //setCursor(*zoomCursor);
        break;
    case CurAction3d_DynamicRotation:
        //setCursor(*rotCursor);
        break;
    case CurAction3d_GlobalPanning:
        setCursor(*globPanCursor);
        break;
    case CurAction3d_WindowZooming:
        setCursor(*handCursor);
        break;
    case CurAction3d_Nothing:
    default:
        setCursor(*defCursor);
        break;
    }

}

void OccWidget::onLButtonDown(const int/*Qt::MouseButtons*/ nFlags, const QPoint point)
{

    if (nFlags & MULTISELECTIONKEY)
        m_context->ShiftSelect();
    else
        m_context->Select();

    emit onSelectionChanged();

    // save the current mouse coordinate in min
    m_xmin = point.x();
    m_ymin = point.y();
    m_xmax = point.x();
    m_ymax = point.y();

    if (m_hlrModeIsOn)
        m_view->SetComputedMode(Standard_False);

    m_currentMode = CurAction3d_DynamicRotation;
    m_view->StartRotation(point.x(), point.y(), 0.4);

    activateCursor(m_currentMode);

}

void OccWidget::onMButtonDown(const int/*Qt::MouseButtons*/ nFlags, const QPoint point)
{

}

void OccWidget::onRButtonDown(const int/*Qt::MouseButtons*/ nFlags, const QPoint point)
{

    // save the current mouse coordinate in min
    m_xmin = point.x();
    m_ymin = point.y();
    m_xmax = point.x();
    m_ymax = point.y();

    m_currentMode = CurAction3d_DynamicPanning;        
    activateCursor(m_currentMode);

}

void OccWidget::onLButtonUp(Qt::MouseButtons nFlags, const QPoint point)
{

    m_currentMode = CurAction3d_Nothing;        
    activateCursor(m_currentMode);

}

void OccWidget::onMButtonUp(Qt::MouseButtons /*nFlags*/, const QPoint /*point*/)
{

    m_currentMode = CurAction3d_Nothing;        
    activateCursor(m_currentMode);

}

void OccWidget::onRButtonUp(Qt::MouseButtons /*nFlags*/, const QPoint point)
{

    m_currentMode = CurAction3d_Nothing;        
    activateCursor(m_currentMode);

}

void OccWidget::onMouseMove(Qt::MouseButtons nFlags, const QPoint point)
{

    if (nFlags & Qt::LeftButton || nFlags & Qt::RightButton || nFlags & Qt::MidButton)
    {
        switch (m_currentMode)
        {
        case CurAction3d_Nothing:
            break;
        case CurAction3d_WindowSelection:
            m_xmax = point.x();
            m_ymax = point.y();
            DrawRectangle(m_xmin, m_ymin, m_xmax, m_ymax, Standard_False);
            DrawRectangle(m_xmin, m_ymin, m_xmax, m_ymax, Standard_True);
            m_context->Select(m_xmin, m_ymin, m_xmax, m_ymax, m_view);
            emit onSelectionChanged();
            break;
        case CurAction3d_DynamicZooming:
            break;
        case CurAction3d_WindowZooming:
            m_xmax = point.x();
            m_ymax = point.y();
            DrawRectangle(m_xmin, m_ymin, m_xmax, m_ymax, Standard_False);
            DrawRectangle(m_xmin, m_ymin, m_xmax, m_ymax, Standard_True);
            break;
        case CurAction3d_DynamicPanning:
            m_view->Pan(point.x() - m_xmax, m_ymax - point.y());
            m_xmax = point.x();
            m_ymax = point.y();
            break;
        case CurAction3d_GlobalPanning:
            break;
        case CurAction3d_DynamicRotation:
            m_view->Rotation(point.x(), point.y());
            m_view->Redraw();
            break;
        default:
            Standard_Failure::Raise("incompatible Current Mode");
            break;
        }
    }

    m_context->MoveTo(point.x(), point.y(), m_view);

}

void OccWidget::Popup(const int /*x*/, const int /*y*/)
{

}

void OccWidget::addItemInPopup(QMenu* /*theMenu*/)
{
}

void OccWidget::DrawRectangle(const int MinX, const int MinY, const int MaxX, const int MaxY, const bool Draw)
{ 

    static Standard_Integer StoredMinX, StoredMaxX, StoredMinY, StoredMaxY;
    static Standard_Boolean m_IsVisible;

    StoredMinX = (MinX < MaxX) ? MinX: MaxX ;
    StoredMinY = (MinY < MaxY) ? MinY: MaxY ;
    StoredMaxX = (MinX > MaxX) ? MinX: MaxX ;
    StoredMaxY = (MinY > MaxY) ? MinY: MaxY ;

    QRect aRect;
    aRect.setRect(StoredMinX, StoredMinY, abs(StoredMaxX-StoredMinX), abs(StoredMaxY-StoredMinY));

    if (!m_rectBand) 
    {
        m_rectBand = new QRubberBand(QRubberBand::Rectangle, this);
        m_rectBand->setStyle(new QWindowsStyle);
        m_rectBand->setGeometry(aRect);
        m_rectBand->show();
    }

    if (m_IsVisible && !Draw) // move or up  : erase at the old position
    {
        m_rectBand->hide();
        delete m_rectBand;
        m_rectBand = 0;
        m_IsVisible = false;
    }

    if (Draw) // move : draw
    {
        //aRect.setRect(StoredMinX, StoredMinY, abs(StoredMaxX-StoredMinX), abs(StoredMaxY-StoredMinY));
        m_IsVisible = true;
        m_rectBand->setGeometry(aRect);
        //m_rectBand->show();
    }

}

void OccWidget::noActiveActions()
{
    for (int i = ViewFitAllId; i < ViewHlrOffId ; ++i)
    {
        QAction* anAction = m_viewActions->at(i);
        if((anAction == m_viewActions->at(ViewFitAreaId)) ||
            (anAction == m_viewActions->at(ViewZoomId)) ||
            (anAction == m_viewActions->at(ViewPanId)) ||
            (anAction == m_viewActions->at(ViewGlobalPanId)) ||
            (anAction == m_viewActions->at(ViewRotationId)))
        {
            setCursor(*defCursor);
            anAction->setCheckable(true);
            anAction->setChecked(false);
        }
    }
}

void OccWidget::onBackground()
{
    QColor aColor ;
    Standard_Real R1;
    Standard_Real G1;
    Standard_Real B1;
    m_view->BackgroundColor(Quantity_TOC_RGB,R1,G1,B1);
    aColor.setRgb(R1*255,G1*255,B1*255);

    QColor aRetColor = QColorDialog::getColor(aColor);

    if(aRetColor.isValid())
    {
        R1 = aRetColor.red()/255.;
        G1 = aRetColor.green()/255.;
        B1 = aRetColor.blue()/255.;
        m_view->SetBackgroundColor(Quantity_TOC_RGB,R1,G1,B1);
    }
    m_view->Redraw();
}

bool OccWidget::dump(Standard_CString theFile)
{

    m_view->Redraw();
    QString ext = QFileInfo(QString(theFile)).completeSuffix();
    if (!ext.compare("ps") || !ext.compare("eps") || !ext.compare("tex") || !ext.compare("pdf") || !ext.compare("svg") || !ext.compare("pgf"))
    {
        Graphic3d_ExportFormat exFormat;
        if (!ext.compare("ps"))
            exFormat = Graphic3d_EF_PostScript;
        if (!ext.compare("eps"))
            exFormat = Graphic3d_EF_EnhPostScript;
        if (!ext.compare("tex"))
            exFormat = Graphic3d_EF_TEX;
        if (!ext.compare("pdf"))
            exFormat = Graphic3d_EF_PDF;
        if (!ext.compare("svg"))
            exFormat = Graphic3d_EF_SVG;
        if (!ext.compare("pgf"))
            exFormat = Graphic3d_EF_PGF;
        try
        {
            m_view->View()->Export(theFile, exFormat);
        }
        catch(...)
        {
            return false;
        }
        return true;
    }
    return m_view->Dump(theFile);

}

Handle(V3d_View)& OccWidget::getView()
{

    return m_view;

}

Handle(AIS_InteractiveContext)& OccWidget::getContext()
{

    return m_context;

}

OccWidget::CurrentAction3d OccWidget::getCurrentMode()
{

    return m_currentMode;

}



