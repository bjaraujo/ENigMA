#pragma once

#include <QtGui/QWidget>
#include <QtGui/QAction>
#include <QtCore/QList>

#include <AIS_InteractiveContext.hxx>
#include <AIS_InteractiveObject.hxx>
#include <AIS_Shape.hxx>

#include <OpenGl_GraphicDriver.hxx>

#include <Graphic3d_NameOfMaterial.hxx>
#include <Graphic3d_GraphicDriver.hxx>

#include <Aspect_DisplayConnection.hxx>
#include <TCollection_AsciiString.hxx>

#include <V3d_Viewer.hxx>
#include <V3d_View.hxx>

class TopoDS_Shape;
class QRubberBand;

//class COMMONSAMPLE_EXPORT View: public QWidget
class OccWidget: public QWidget
{
    Q_OBJECT
protected:
    enum CurrentAction3d {CurAction3d_Nothing, CurAction3d_WindowSelection, 
                          CurAction3d_DynamicZooming, CurAction3d_WindowZooming, 
                          CurAction3d_DynamicPanning, CurAction3d_GlobalPanning, CurAction3d_DynamicRotation};

public:
    enum ViewAction {ViewFitAllId, ViewFitAreaId, ViewZoomId, ViewPanId, ViewGlobalPanId,
                     ViewFrontId, ViewBackId, ViewTopId, ViewBottomId, ViewLeftId, ViewRightId,
                     ViewAxoId, ViewRotationId, ViewResetId, ViewHlrOffId, ViewHlrOnId};

    OccWidget();
    ~OccWidget();

    virtual void                  setContext(Handle(AIS_InteractiveContext) theContext);
    virtual void                  init();
    bool                          dump(Standard_CString theFile);
    QList<QAction*>*              getViewActions();
    void                          noActiveActions();
    bool                          isShadingMode();

    static QString                getMessages(int type, TopAbs_ShapeEnum aSubShapeType, TopAbs_ShapeEnum aShapeType);
    static QString                getShapeType(TopAbs_ShapeEnum aShapeType);

    Standard_EXPORT static void   onButtonUserAction(int ExerciceSTEP, Handle(AIS_InteractiveContext)&);
    Standard_EXPORT static void   doSelection(int Id, Handle(AIS_InteractiveContext)&);
    Standard_EXPORT static void   onSetSelectionMode(Handle(AIS_InteractiveContext)&, Standard_Integer&, TopAbs_ShapeEnum& SelectionMode, Standard_Boolean&);

    virtual QPaintEngine*         paintEngine() const;

signals:
    void                          onSelectionChanged();

public slots:
    void fitAll();
    void fitArea();
    void zoom();
    void pan();
    void globalPan();
    void front();
    void back();
    void top();
    void bottom();
    void left();
    void right();
    void axo();
    void rotation();
    void reset();
    void hlrOn();
    void hlrOff();
    void updateToggled(bool);
    void onBackground();

protected:
    virtual void paintEvent(QPaintEvent*);
    virtual void resizeEvent(QResizeEvent*);
    virtual void mousePressEvent(QMouseEvent*);
    virtual void mouseReleaseEvent(QMouseEvent*);
    virtual void mouseMoveEvent(QMouseEvent*);
    virtual void wheelEvent(QWheelEvent*);

    virtual void                  addItemInPopup(QMenu*);

    Handle(V3d_View)&                     getView();
    Handle(AIS_InteractiveContext)&       getContext();
    void                                  activateCursor(const CurrentAction3d);
    void                                  Popup(const int x, const int y);
    CurrentAction3d                          getCurrentMode();

    virtual void onLButtonDown(const int nFlags, const QPoint point);
    virtual void onMButtonDown(const int nFlags, const QPoint point);
    virtual void onRButtonDown(const int nFlags, const QPoint point);
    virtual void onLButtonUp(Qt::MouseButtons nFlags, const QPoint point);
    virtual void onMButtonUp(Qt::MouseButtons nFlags, const QPoint point);
    virtual void onRButtonUp(Qt::MouseButtons nFlags, const QPoint point);
    virtual void onMouseMove(Qt::MouseButtons nFlags, const QPoint point);

private:
    void initCursors();
    void initViewActions();
    void DrawRectangle(const int MinX, const int MinY, const int MaxX, const int MaxY, const bool Draw);

private:
    bool                            m_first;
    bool                            m_drawRect;           // set when a rect is used for selection or magnify 
    Handle(V3d_View)                m_view;
    Handle(AIS_InteractiveContext)  m_context;
    CurrentAction3d                 m_currentMode;
    Standard_Integer                m_xmin;
    Standard_Integer                m_ymin;
    Standard_Integer                m_xmax;
    Standard_Integer                m_ymax;
    Quantity_Factor                 m_curZoom;
    Standard_Boolean                m_hlrModeIsOn;
    QList<QAction*>*                m_viewActions;
    QRubberBand*                    m_rectBand; //!< selection rectangle rubber band
};


