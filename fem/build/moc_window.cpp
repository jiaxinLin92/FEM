/****************************************************************************
** Meta object code from reading C++ file 'window.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../window.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'window.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_Window_t {
    QByteArrayData data[17];
    char stringdata[234];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Window_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Window_t qt_meta_stringdata_Window = {
    {
QT_MOC_LITERAL(0, 0, 6), // "Window"
QT_MOC_LITERAL(1, 7, 22), // "signGuiViewportChanged"
QT_MOC_LITERAL(2, 30, 0), // ""
QT_MOC_LITERAL(3, 31, 4), // "name"
QT_MOC_LITERAL(4, 36, 12), // "new_geometry"
QT_MOC_LITERAL(5, 49, 14), // "signFrameReady"
QT_MOC_LITERAL(6, 64, 16), // "signMousePressed"
QT_MOC_LITERAL(7, 81, 12), // "QMouseEvent*"
QT_MOC_LITERAL(8, 94, 5), // "event"
QT_MOC_LITERAL(9, 100, 17), // "signMouseReleased"
QT_MOC_LITERAL(10, 118, 22), // "signMouseDoubleClicked"
QT_MOC_LITERAL(11, 141, 14), // "signMouseMoved"
QT_MOC_LITERAL(12, 156, 14), // "signKeyPressed"
QT_MOC_LITERAL(13, 171, 10), // "QKeyEvent*"
QT_MOC_LITERAL(14, 182, 15), // "signKeyReleased"
QT_MOC_LITERAL(15, 198, 22), // "signWheelEventOccurred"
QT_MOC_LITERAL(16, 221, 12) // "QWheelEvent*"

    },
    "Window\0signGuiViewportChanged\0\0name\0"
    "new_geometry\0signFrameReady\0"
    "signMousePressed\0QMouseEvent*\0event\0"
    "signMouseReleased\0signMouseDoubleClicked\0"
    "signMouseMoved\0signKeyPressed\0QKeyEvent*\0"
    "signKeyReleased\0signWheelEventOccurred\0"
    "QWheelEvent*"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Window[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       9,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    2,   59,    2, 0x06 /* Public */,
       5,    0,   64,    2, 0x06 /* Public */,
       6,    2,   65,    2, 0x06 /* Public */,
       9,    2,   70,    2, 0x06 /* Public */,
      10,    2,   75,    2, 0x06 /* Public */,
      11,    2,   80,    2, 0x06 /* Public */,
      12,    2,   85,    2, 0x06 /* Public */,
      14,    2,   90,    2, 0x06 /* Public */,
      15,    2,   95,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString, QMetaType::QRectF,    3,    4,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    3,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    3,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    3,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 7,    3,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 13,    3,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 13,    3,    8,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 16,    3,    8,

       0        // eod
};

void Window::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Window *_t = static_cast<Window *>(_o);
        switch (_id) {
        case 0: _t->signGuiViewportChanged((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< const QRectF(*)>(_a[2]))); break;
        case 1: _t->signFrameReady(); break;
        case 2: _t->signMousePressed((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 3: _t->signMouseReleased((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 4: _t->signMouseDoubleClicked((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 5: _t->signMouseMoved((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 6: _t->signKeyPressed((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QKeyEvent*(*)>(_a[2]))); break;
        case 7: _t->signKeyReleased((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QKeyEvent*(*)>(_a[2]))); break;
        case 8: _t->signWheelEventOccurred((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QWheelEvent*(*)>(_a[2]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (Window::*_t)(const QString & , const QRectF & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signGuiViewportChanged)) {
                *result = 0;
            }
        }
        {
            typedef void (Window::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signFrameReady)) {
                *result = 1;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signMousePressed)) {
                *result = 2;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signMouseReleased)) {
                *result = 3;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signMouseDoubleClicked)) {
                *result = 4;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signMouseMoved)) {
                *result = 5;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QKeyEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signKeyPressed)) {
                *result = 6;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QKeyEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signKeyReleased)) {
                *result = 7;
            }
        }
        {
            typedef void (Window::*_t)(const QString & , QWheelEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&Window::signWheelEventOccurred)) {
                *result = 8;
            }
        }
    }
}

const QMetaObject Window::staticMetaObject = {
    { &QQuickView::staticMetaObject, qt_meta_stringdata_Window.data,
      qt_meta_data_Window,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *Window::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Window::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_Window.stringdata))
        return static_cast<void*>(const_cast< Window*>(this));
    return QQuickView::qt_metacast(_clname);
}

int Window::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QQuickView::qt_metacall(_c, _id, _a);
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
void Window::signGuiViewportChanged(const QString & _t1, const QRectF & _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void Window::signFrameReady()
{
    QMetaObject::activate(this, &staticMetaObject, 1, Q_NULLPTR);
}

// SIGNAL 2
void Window::signMousePressed(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void Window::signMouseReleased(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void Window::signMouseDoubleClicked(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void Window::signMouseMoved(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void Window::signKeyPressed(const QString & _t1, QKeyEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void Window::signKeyReleased(const QString & _t1, QKeyEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}

// SIGNAL 8
void Window::signWheelEventOccurred(const QString & _t1, QWheelEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}
QT_END_MOC_NAMESPACE
