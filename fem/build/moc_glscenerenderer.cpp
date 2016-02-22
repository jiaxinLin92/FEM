/****************************************************************************
** Meta object code from reading C++ file 'glscenerenderer.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../glscenerenderer.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'glscenerenderer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_Private__Renderer_t {
    QByteArrayData data[9];
    char stringdata[50];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Private__Renderer_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Private__Renderer_t qt_meta_stringdata_Private__Renderer = {
    {
QT_MOC_LITERAL(0, 0, 17), // "Private::Renderer"
QT_MOC_LITERAL(1, 18, 5), // "paint"
QT_MOC_LITERAL(2, 24, 0), // ""
QT_MOC_LITERAL(3, 25, 11), // "setViewport"
QT_MOC_LITERAL(4, 37, 4), // "rect"
QT_MOC_LITERAL(5, 42, 1), // "x"
QT_MOC_LITERAL(6, 44, 1), // "y"
QT_MOC_LITERAL(7, 46, 1), // "w"
QT_MOC_LITERAL(8, 48, 1) // "h"

    },
    "Private::Renderer\0paint\0\0setViewport\0"
    "rect\0x\0y\0w\0h"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Private__Renderer[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   29,    2, 0x0a /* Public */,
       3,    1,   30,    2, 0x0a /* Public */,
       3,    4,   33,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::QRectF,    4,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, QMetaType::Int, QMetaType::Int,    5,    6,    7,    8,

       0        // eod
};

void Private::Renderer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Renderer *_t = static_cast<Renderer *>(_o);
        switch (_id) {
        case 0: _t->paint(); break;
        case 1: _t->setViewport((*reinterpret_cast< const QRectF(*)>(_a[1]))); break;
        case 2: _t->setViewport((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4]))); break;
        default: ;
        }
    }
}

const QMetaObject Private::Renderer::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_Private__Renderer.data,
      qt_meta_data_Private__Renderer,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *Private::Renderer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Private::Renderer::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_Private__Renderer.stringdata))
        return static_cast<void*>(const_cast< Renderer*>(this));
    return QObject::qt_metacast(_clname);
}

int Private::Renderer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}
struct qt_meta_stringdata_GLSceneRenderer_t {
    QByteArrayData data[37];
    char stringdata[486];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GLSceneRenderer_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GLSceneRenderer_t qt_meta_stringdata_GLSceneRenderer = {
    {
QT_MOC_LITERAL(0, 0, 15), // "GLSceneRenderer"
QT_MOC_LITERAL(1, 16, 19), // "signViewportChanged"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 4), // "name"
QT_MOC_LITERAL(4, 42, 4), // "size"
QT_MOC_LITERAL(5, 47, 16), // "signMousePressed"
QT_MOC_LITERAL(6, 64, 12), // "QMouseEvent*"
QT_MOC_LITERAL(7, 77, 5), // "event"
QT_MOC_LITERAL(8, 83, 17), // "signMouseReleased"
QT_MOC_LITERAL(9, 101, 22), // "signMouseDoubleClicked"
QT_MOC_LITERAL(10, 124, 14), // "signMouseMoved"
QT_MOC_LITERAL(11, 139, 14), // "signKeyPressed"
QT_MOC_LITERAL(12, 154, 10), // "QKeyEvent*"
QT_MOC_LITERAL(13, 165, 15), // "signKeyReleased"
QT_MOC_LITERAL(14, 181, 22), // "signWheelEventOccurred"
QT_MOC_LITERAL(15, 204, 12), // "QWheelEvent*"
QT_MOC_LITERAL(16, 217, 4), // "sync"
QT_MOC_LITERAL(17, 222, 7), // "cleanup"
QT_MOC_LITERAL(18, 230, 10), // "itemChange"
QT_MOC_LITERAL(19, 241, 10), // "ItemChange"
QT_MOC_LITERAL(20, 252, 6), // "change"
QT_MOC_LITERAL(21, 259, 14), // "ItemChangeData"
QT_MOC_LITERAL(22, 274, 5), // "value"
QT_MOC_LITERAL(23, 280, 15), // "updatePaintNode"
QT_MOC_LITERAL(24, 296, 8), // "QSGNode*"
QT_MOC_LITERAL(25, 305, 20), // "UpdatePaintNodeData*"
QT_MOC_LITERAL(26, 326, 15), // "mousePressEvent"
QT_MOC_LITERAL(27, 342, 17), // "mouseReleaseEvent"
QT_MOC_LITERAL(28, 360, 21), // "mouseDoubleClickEvent"
QT_MOC_LITERAL(29, 382, 14), // "mouseMoveEvent"
QT_MOC_LITERAL(30, 397, 13), // "keyPressEvent"
QT_MOC_LITERAL(31, 411, 15), // "keyReleaseEvent"
QT_MOC_LITERAL(32, 427, 10), // "wheelEvent"
QT_MOC_LITERAL(33, 438, 19), // "handleWindowChanged"
QT_MOC_LITERAL(34, 458, 13), // "QQuickWindow*"
QT_MOC_LITERAL(35, 472, 6), // "window"
QT_MOC_LITERAL(36, 479, 6) // "paused"

    },
    "GLSceneRenderer\0signViewportChanged\0"
    "\0name\0size\0signMousePressed\0QMouseEvent*\0"
    "event\0signMouseReleased\0signMouseDoubleClicked\0"
    "signMouseMoved\0signKeyPressed\0QKeyEvent*\0"
    "signKeyReleased\0signWheelEventOccurred\0"
    "QWheelEvent*\0sync\0cleanup\0itemChange\0"
    "ItemChange\0change\0ItemChangeData\0value\0"
    "updatePaintNode\0QSGNode*\0UpdatePaintNodeData*\0"
    "mousePressEvent\0mouseReleaseEvent\0"
    "mouseDoubleClickEvent\0mouseMoveEvent\0"
    "keyPressEvent\0keyReleaseEvent\0wheelEvent\0"
    "handleWindowChanged\0QQuickWindow*\0"
    "window\0paused"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GLSceneRenderer[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      20,   14, // methods
       2,  190, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       8,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    2,  114,    2, 0x06 /* Public */,
       5,    2,  119,    2, 0x06 /* Public */,
       8,    2,  124,    2, 0x06 /* Public */,
       9,    2,  129,    2, 0x06 /* Public */,
      10,    2,  134,    2, 0x06 /* Public */,
      11,    2,  139,    2, 0x06 /* Public */,
      13,    2,  144,    2, 0x06 /* Public */,
      14,    2,  149,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      16,    0,  154,    2, 0x0a /* Public */,
      17,    0,  155,    2, 0x0a /* Public */,
      18,    2,  156,    2, 0x09 /* Protected */,
      23,    2,  161,    2, 0x09 /* Protected */,
      26,    1,  166,    2, 0x09 /* Protected */,
      27,    1,  169,    2, 0x09 /* Protected */,
      28,    1,  172,    2, 0x09 /* Protected */,
      29,    1,  175,    2, 0x09 /* Protected */,
      30,    1,  178,    2, 0x09 /* Protected */,
      31,    1,  181,    2, 0x09 /* Protected */,
      32,    1,  184,    2, 0x09 /* Protected */,
      33,    1,  187,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString, QMetaType::QRectF,    3,    4,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 6,    3,    7,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 6,    3,    7,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 6,    3,    7,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 6,    3,    7,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 12,    3,    7,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 12,    3,    7,
    QMetaType::Void, QMetaType::QString, 0x80000000 | 15,    3,    7,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 19, 0x80000000 | 21,   20,   22,
    0x80000000 | 24, 0x80000000 | 24, 0x80000000 | 25,    2,    2,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, 0x80000000 | 12,    7,
    QMetaType::Void, 0x80000000 | 12,    7,
    QMetaType::Void, 0x80000000 | 15,    7,
    QMetaType::Void, 0x80000000 | 34,   35,

 // properties: name, type, flags
       3, QMetaType::QString, 0x00095003,
      36, QMetaType::Bool, 0x00095103,

       0        // eod
};

void GLSceneRenderer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GLSceneRenderer *_t = static_cast<GLSceneRenderer *>(_o);
        switch (_id) {
        case 0: _t->signViewportChanged((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< const QRectF(*)>(_a[2]))); break;
        case 1: _t->signMousePressed((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 2: _t->signMouseReleased((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 3: _t->signMouseDoubleClicked((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 4: _t->signMouseMoved((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QMouseEvent*(*)>(_a[2]))); break;
        case 5: _t->signKeyPressed((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QKeyEvent*(*)>(_a[2]))); break;
        case 6: _t->signKeyReleased((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QKeyEvent*(*)>(_a[2]))); break;
        case 7: _t->signWheelEventOccurred((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< QWheelEvent*(*)>(_a[2]))); break;
        case 8: _t->sync(); break;
        case 9: _t->cleanup(); break;
        case 10: _t->itemChange((*reinterpret_cast< ItemChange(*)>(_a[1])),(*reinterpret_cast< const ItemChangeData(*)>(_a[2]))); break;
        case 11: { QSGNode* _r = _t->updatePaintNode((*reinterpret_cast< QSGNode*(*)>(_a[1])),(*reinterpret_cast< UpdatePaintNodeData*(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< QSGNode**>(_a[0]) = _r; }  break;
        case 12: _t->mousePressEvent((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 13: _t->mouseReleaseEvent((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 14: _t->mouseDoubleClickEvent((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 15: _t->mouseMoveEvent((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 16: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 17: _t->keyReleaseEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 18: _t->wheelEvent((*reinterpret_cast< QWheelEvent*(*)>(_a[1]))); break;
        case 19: _t->handleWindowChanged((*reinterpret_cast< QQuickWindow*(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , const QRectF & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signViewportChanged)) {
                *result = 0;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signMousePressed)) {
                *result = 1;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signMouseReleased)) {
                *result = 2;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signMouseDoubleClicked)) {
                *result = 3;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QMouseEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signMouseMoved)) {
                *result = 4;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QKeyEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signKeyPressed)) {
                *result = 5;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QKeyEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signKeyReleased)) {
                *result = 6;
            }
        }
        {
            typedef void (GLSceneRenderer::*_t)(const QString & , QWheelEvent * );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLSceneRenderer::signWheelEventOccurred)) {
                *result = 7;
            }
        }
    }
}

const QMetaObject GLSceneRenderer::staticMetaObject = {
    { &QQuickItem::staticMetaObject, qt_meta_stringdata_GLSceneRenderer.data,
      qt_meta_data_GLSceneRenderer,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *GLSceneRenderer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GLSceneRenderer::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_GLSceneRenderer.stringdata))
        return static_cast<void*>(const_cast< GLSceneRenderer*>(this));
    return QQuickItem::qt_metacast(_clname);
}

int GLSceneRenderer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QQuickItem::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 20)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 20;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 20)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 20;
    }
#ifndef QT_NO_PROPERTIES
      else if (_c == QMetaObject::ReadProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< QString*>(_v) = getTexName(); break;
        case 1: *reinterpret_cast< bool*>(_v) = isPaused(); break;
        default: break;
        }
        _id -= 2;
    } else if (_c == QMetaObject::WriteProperty) {
        void *_v = _a[0];
        switch (_id) {
        case 0: setTexName(*reinterpret_cast< QString*>(_v)); break;
        case 1: setPaused(*reinterpret_cast< bool*>(_v)); break;
        default: break;
        }
        _id -= 2;
    } else if (_c == QMetaObject::ResetProperty) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 2;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 2;
    } else if (_c == QMetaObject::RegisterPropertyMetaType) {
        if (_id < 2)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 2;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void GLSceneRenderer::signViewportChanged(const QString & _t1, const QRectF & _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void GLSceneRenderer::signMousePressed(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void GLSceneRenderer::signMouseReleased(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void GLSceneRenderer::signMouseDoubleClicked(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void GLSceneRenderer::signMouseMoved(const QString & _t1, QMouseEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void GLSceneRenderer::signKeyPressed(const QString & _t1, QKeyEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void GLSceneRenderer::signKeyReleased(const QString & _t1, QKeyEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void GLSceneRenderer::signWheelEventOccurred(const QString & _t1, QWheelEvent * _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}
QT_END_MOC_NAMESPACE
