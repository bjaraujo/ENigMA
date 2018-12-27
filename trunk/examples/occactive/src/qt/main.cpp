
// QT includes
#include <QtGui/QApplication>
#include <QtGui/QCleanlooksStyle>

#include "occ/OccManager.h"

#include "MainWindow.h"

QSharedPointer<OccManager> OccManager::m_instance;

int main(int argc, char** argv)
{
  
    try 
    {

        // QT Stuff
        QApplication app(argc, argv);

        QtMainWindow mMainWindow(&app);

        QSharedPointer<OccManager> aManager = OccManager::instance();

        mMainWindow.setManager(aManager);

        mMainWindow.show();

        return app.exec();

    }
    catch(...)
    {
        return -1;
    }

}
