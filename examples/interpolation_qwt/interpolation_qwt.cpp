#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QPushButton>
#include <qwt/qwt.h>
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_point_data.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_zoomer.h>
#include <qwt/qwt_plot_panner.h>
#include <QPen>
#include <QVector>
#include <algorithm>
#include <cmath>

#include <ALFI.h>

template <typename T>
QVector<T> to_qvector(const std::vector<T>& v) {
    return QVector<T>(v.begin(), v.end());
}

double f(double x) {
    return -3 * std::sin(10 * x) + 10 * sin(std::abs(x) + 0.5 * x);
}

template <typename Container>
Container f(const Container& X) {
    Container Y;
    Y.reserve(X.size());
    for (const auto& x : X) {
        Y.push_back(f(x));
    }
    return Y;
}

class PlotWindow final : public QWidget {
public:
    PlotWindow() {
        // Инициализация интерфейса
        static const QStringList distribution_types {
            "Uniform",
            "Chebyshev", "Chebyshev Stretched", "Chebyshev Ellipse", "Chebyshev Ellipse Stretched",
            "Circle Projection", "Ellipse Projection",
            "Sigmoid", "Stretched Sigmoid",
            "Error Function", "Stretched Error Function"
        };

        _plot = new QwtPlot();
        _plot->setCanvasBackground(Qt::white);
        QwtPlotGrid* grid = new QwtPlotGrid();
        grid->attach(_plot);

        QwtPlotZoomer* zoomer = new QwtPlotZoomer(QwtPlot::xBottom, QwtPlot::yLeft, _plot->canvas());
        zoomer->setRubberBand(QwtPlotZoomer::NoRubberBand);
        zoomer->setTrackerMode(QwtPlotZoomer::AlwaysOff);
        zoomer->setZoomBase(true);
        zoomer->setMousePattern(QwtEventPattern::MouseSelect2, Qt::MidButton);

        QwtPlotPanner* panner = new QwtPlotPanner(_plot->canvas());
        panner->setMouseButton(Qt::LeftButton);
        
        _control_panel = new QWidget();
        QVBoxLayout* control_layout = new QVBoxLayout(_control_panel);

        _n_spin_box = create_spin_box<int, QSpinBox>(11, 1);
        _nn_spin_box = create_spin_box<int, QSpinBox>(10000, 1);
        _a_spin_box = create_spin_box<double, QDoubleSpinBox>(-10, -0.1);
        _b_spin_box = create_spin_box<double, QDoubleSpinBox>(10, 0.1);

        _distribution_combo = new QComboBox();
        _distribution_combo->addItems(distribution_types);

        _function_checkbox = new QCheckBox("Function");
        _points_checkbox = new QCheckBox("Points");
        _lagrange_checkbox = new QCheckBox("Lagrange Polynomial");
        _newton_checkbox = new QCheckBox("Newton Polynomial");
        _barycentric_checkbox = new QCheckBox("Barycentric Formula");
        _step_spline_checkbox = new QCheckBox("Step Spline");
        _linear_spline_checkbox = new QCheckBox("Linear Spline");
        _quadratic_spline_checkbox = new QCheckBox("Quadratic Spline");

        _function_checkbox->setChecked(true);
        _points_checkbox->setChecked(true);
        _barycentric_checkbox->setChecked(true);

        _view_reset_button = new QPushButton("Reset View");

        // Подключаем сигналы и слоты
        connect(_n_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);
        connect(_nn_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);

        connect(_a_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot);
        connect(_b_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot);

        for (const auto checkbox : {_function_checkbox, _points_checkbox, _lagrange_checkbox,
            _newton_checkbox, _barycentric_checkbox, _step_spline_checkbox,
            _linear_spline_checkbox, _quadratic_spline_checkbox}) {
            connect(checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);
        }

        for (const auto combo : {_distribution_combo}) {
            connect(combo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &PlotWindow::update_plot);
        }

        connect(_view_reset_button, &QPushButton::clicked, this, [this]() {
            _default_axis_ranges = true;
            update_plot();
        });

        // Создание элементов управления
        QVBoxLayout* n_layout = new QVBoxLayout();
        QHBoxLayout* n_label_spin_box_layout = new QHBoxLayout();
        n_label_spin_box_layout->addWidget(new QLabel("N:"));
        n_label_spin_box_layout->addWidget(_n_spin_box);
        n_layout->addLayout(n_label_spin_box_layout);

        QHBoxLayout* nn_label_spin_box_layout = new QHBoxLayout();
        nn_label_spin_box_layout->addWidget(new QLabel("nn:"));
        nn_label_spin_box_layout->addWidget(_nn_spin_box);
        n_layout->addLayout(nn_label_spin_box_layout);

        QHBoxLayout* a_b_layout = new QHBoxLayout();
        QLabel* a_label = new QLabel("a:");
        QLabel* b_label = new QLabel("b:");
        a_b_layout->addWidget(a_label);
        a_b_layout->addWidget(_a_spin_box);
        a_b_layout->addWidget(b_label);
        a_b_layout->addWidget(_b_spin_box);

        QVBoxLayout* interpolation_layout = new QVBoxLayout();
        interpolation_layout->addWidget(_function_checkbox);
        interpolation_layout->addWidget(_points_checkbox);
        interpolation_layout->addWidget(_lagrange_checkbox);
        interpolation_layout->addWidget(_newton_checkbox);
        interpolation_layout->addWidget(_barycentric_checkbox);
        interpolation_layout->addWidget(_step_spline_checkbox);
        interpolation_layout->addWidget(_linear_spline_checkbox);
        interpolation_layout->addWidget(_quadratic_spline_checkbox);

        control_layout->addLayout(n_layout);
        control_layout->addLayout(a_b_layout);
        control_layout->addLayout(interpolation_layout);
        control_layout->addWidget(_distribution_combo);
        control_layout->addWidget(_view_reset_button);

        _control_panel->setMinimumWidth(250);

        QHBoxLayout* main_layout = new QHBoxLayout();
        main_layout->addWidget(_plot, 1);
        main_layout->addWidget(_control_panel);

        setLayout(main_layout);

        update_plot();
    }

private:
    template <typename Number, typename SpinBoxType>
    static SpinBoxType* create_spin_box(Number initial, Number step) {
        SpinBoxType* spin_box = new SpinBoxType();
        spin_box->setRange(std::min(std::numeric_limits<Number>::lowest(), -std::numeric_limits<Number>::infinity()),
                           std::max(std::numeric_limits<Number>::max(), std::numeric_limits<Number>::infinity()));
        spin_box->setValue(initial);
        spin_box->setSingleStep(step);
        return spin_box;
    }

    void update_plot() {
        static const std::vector<QColor> colors = {
            QColor(0, 102, 204), QColor(204, 0, 0), QColor(0, 153, 0),
            QColor(204, 102, 0), QColor(102, 0, 153), QColor(30, 30, 30),
            QColor(0, 128, 128), QColor(204, 0, 102)
        };

        const size_t N = _n_spin_box->value();
        const size_t nn = _nn_spin_box->value();
        const double a = _a_spin_box->value();
        const double b = _b_spin_box->value();
        const int dist_type = _distribution_combo->currentIndex() + 1;

        std::vector<double> X;
        switch (static_cast<alfi::dist::Type>(dist_type)) {
            case alfi::dist::Type::UNIFORM: X = alfi::dist::uniform(N, a, b); break;
            case alfi::dist::Type::CHEBYSHEV: X = alfi::dist::chebyshev(N, a, b); break;
            case alfi::dist::Type::CHEBYSHEV_STRETCHED: X = alfi::dist::chebyshev_stretched(N, a, b); break;
            case alfi::dist::Type::CHEBYSHEV_ELLIPSE: X = alfi::dist::chebyshev_ellipse(N, a, b, 2.0); break;
            case alfi::dist::Type::CHEBYSHEV_ELLIPSE_STRETCHED: X = alfi::dist::chebyshev_ellipse_stretched(N, a, b, 2.0); break;
            case alfi::dist::Type::CIRCLE_PROJECTION: X = alfi::dist::circle_proj(N, a, b); break;
            case alfi::dist::Type::ELLIPSE_PROJECTION: X = alfi::dist::ellipse_proj(N, a, b, 2.0); break;
            case alfi::dist::Type::SIGMOID: X = alfi::dist::sigmoid(N, a, b, 16.0); break;
            case alfi::dist::Type::SIGMOID_STRETCHED: X = alfi::dist::sigmoid_stretched(N, a, b, 16.0); break;
            case alfi::dist::Type::ERF: X = alfi::dist::erf(N, a, b, 8.0); break;
            case alfi::dist::Type::ERF_STRETCHED: X = alfi::dist::erf_stretched(N, a, b, 8.0); break;
            default: X = alfi::dist::uniform(N, a, b); break;
        }

        const std::vector<double> Y = f(X);

        _plot->detachItems(QwtPlotItem::Rtti_PlotCurve);

        if (_points_checkbox->isChecked()) {
            QwtPlotCurve* curve = new QwtPlotCurve();
            curve->setStyle(QwtPlotCurve::Dots);
            curve->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QBrush(Qt::red), QPen(Qt::red), QSize(5, 5)));
            curve->setData(new QwtPointArrayData(to_qvector(X), to_qvector(Y)));
            curve->attach(_plot);
        }

        if (_function_checkbox->isChecked()) {
            QwtPlotCurve* curve = new QwtPlotCurve("Function");
            curve->setPen(QPen(Qt::blue));
            curve->setData(new QwtPointArrayData(to_qvector(X), to_qvector(Y)));
            curve->attach(_plot);
        }

        _plot->replot();
    }

private:
    QwtPlot* _plot;
    QWidget* _control_panel;
    QSpinBox* _n_spin_box;
    QSpinBox* _nn_spin_box;
    QDoubleSpinBox* _a_spin_box;
    QDoubleSpinBox* _b_spin_box;
    QComboBox* _distribution_combo;
    QCheckBox* _function_checkbox;
    QCheckBox* _points_checkbox;
    QCheckBox* _lagrange_checkbox;
    QCheckBox* _newton_checkbox;
    QCheckBox* _barycentric_checkbox;
    QCheckBox* _step_spline_checkbox;
    QCheckBox* _linear_spline_checkbox;
    QCheckBox* _quadratic_spline_checkbox;
    QPushButton* _view_reset_button;
    bool _default_axis_ranges = false;
};

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    PlotWindow window;
    window.setWindowTitle("Interpolation");
    window.resize(1280, 720);
    window.show();
    return QApplication::exec();
}