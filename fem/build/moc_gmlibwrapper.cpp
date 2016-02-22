/****************************************************************************
** Meta object code from reading C++ file 'gmlibwrapper.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../gmlibwrapper.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'gmlibwrapper.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_GMlibWrapper_t {
    QByteArrayData data[17];
    char stringdata[210];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GMlibWrapper_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GMlibWrapper_t qt_meta_stringdata_GMlibWrapper = {
    {
QT_MOC_LITERAL(0, 0, 12), // "GMlibWrapper"
QT_MOC_LITERAL(1, 13, 14), // "signFrameReady"
QT_MOC_LITERAL(2, 28, 0), // ""
QT_MOC_LITERAL(3, 29, 20), // "changeRenderGeometry"
QT_MOC_LITERAL(4, 50, 4), // "name"
QT_MOC_LITERAL(5, 55, 12), // "new_geometry"
QT_MOC_LITERAL(6, 68, 12), // "mousePressed"
QT_MOC_LITERAL(7, 81, 12), // "QMouseEvent*"
QT_MOC_LITERAL(8, 94, 5), // "event"
QT_MOC_LITERAL(9, 100, 13), // "mouseReleased"
QT_MOC_LITERAL(10, 114, 18), // "mouseDoubleClicked"
QT_MOC_LITERAL(11, 133, 10), // "mouseMoved"
QT_MOC_LITERAL(12, 144, 10), // "keyPressed"
QT_MOC_LITERAL(13, 155, 10), // "QKeyEvent*"
QT_MOC_LITERAL(14, 166, 11), // "keyReleased"
QT_MOC_LITERAL(15, 178, 18), // "wheelEventOccurred"
QT_MOC_LITERAL(16, 197, 12) // "QWheelEvent*"

    },
    "GMlibWrapper\0signFrameReady\0\0"
    "changeRenderGeometry\0name\0new_geometry\0"
    "mousePressed\0QMouseEvent*\0event\0"
    "mouseReleased\0mouseDoubleClicked\0"
    "mouseMoved\0keyPressed\0QKeyEvent*\0"
    "keyReleased\0wheelEventOccurred\0"
    "QWheelEvent*"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GMlibWrapper[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   59,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    2,   60,    2, 0x0a /* Public */,
       6,    2,   65,    2, 0x0a /* Public */,
       9,    2,   70,    2, 0x0a /* Public */,
      10,    2,   75,    2, 0x0a /* Public */,
      11,    2,   80,    2, 0x0a /* Public */,
      12,    2,   85,    2, 0x0a /* Public */,
      14,    2,   90,    2, 0x0a /* Public */,
      15,    2,   95,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::QString, QMetaType::QRectF,    4,    5,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    4,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    4,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    4,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    4,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 13,    4,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 13,    4,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 16,    4,    8,

       0        // eod
};

void GMlibWrapper::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GMlibWrapper *_t = static_cast<GMlibWrapper *>(_o);
        switch (_id) {
        case 0: _t->signFrameReady(); break;
        case 1: _t->changeRenderGeometry((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< const QRectF(*)>(_a[2]))); break;
        case 2: _t->mousePressed((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 3: _t->mouseReleased((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 4: _t->mouseDoubleClicked((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 5: _t->mouseMoved((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 6: _t->keyPressed((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QKeyEvent*(*)>(_a[2]))); break;
        case 7: _t->keyReleased((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QKeyEvent*(*)>(_a[2]))); break;
        case 8: _t->wheelEventOccurred((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QWheelEvent*(*)>(_a[2]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (GMlibWrapper::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GMlibWrapper::signFrameReady)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject GMlibWrapper::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_GMlibWrapper.data,
      qt_meta_data_GMlibWrapper,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *GMlibWrapper::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GMlibWrapper::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_GMlibWrapper.stringdata))
        return static_cast<void*>(const_cast< GMlibWrapper*>(this));
    return QObject::qt_metacast(_clname);
}

int GMlibWrapper::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}

// SIGNAL 0
void GMlibWrapper::signFrameReady()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
