#include <qcustomplot.h>
#include <ALFI/ratf.h>

#include <iomanip>
#include <iostream>
#include <vector>

constexpr double factorial(int n) {
	double result = 1;
	for (int i = 2; i <= n; ++i) {
		result *= i;
	}
	return result;
}

std::vector<double> maclaurin_exp(int order) {
	std::vector<double> coeffs(order + 1);
	for (int n = 0; n <= order; ++n) {
		coeffs[order-n] = 1.0 / factorial(n);
	}
	return coeffs;
}

std::vector<double> maclaurin_sin(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 0; n <= order; ++n) {
		if (n % 2 == 1) {
			coeffs[order-n] = (n % 4 == 1 ? 1.0 : -1.0) / factorial(n);
		}
	}
	return coeffs;
}

std::vector<double> maclaurin_cos(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 0; n <= order; ++n) {
		if (n % 2 == 0) {
			coeffs[order-n] = (n % 4 == 0 ? 1.0 : -1.0) / factorial(n);
		}
	}
	return coeffs;
}

std::vector<double> maclaurin_ln(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 1; n <= order; ++n) {
		coeffs[order-n] = (n % 2 == 1 ? 1.0 : -1.0) / n;
	}
	return coeffs;
}

std::vector<double> maclaurin_arctan(int order) {
	std::vector<double> coeffs(order + 1, 0);
	for (int n = 0; n <= order; ++n) {
		if (n % 2 == 1) {
			coeffs[order-n] = (n % 4 == 1 ? 1.0 : -1.0) / n;
		}
	}
	return coeffs;
}

class PlotWindow final : public QWidget {
private:
	template <typename Number>
	struct VariableParams {
		Number initial;
		Number min;
		Number max;
		Number step;
	};

public:
	PlotWindow() {
		static const VariableParams<int> order_params = {10, 1, 40, 1};
		static const VariableParams<int> n_params = {4, 1, 10, 1};
		static const VariableParams<int> m_params = {4, 1, 10, 1};

		_plot = new QCustomPlot();
		_plot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
		_plot->legend->setVisible(true);
		_plot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));

		_control_panel = new QWidget();
		QVBoxLayout* control_layout = new QVBoxLayout(_control_panel);

		_order_spin_box = create_spin_box<int, QSpinBox>(order_params);
		_n_spin_box = create_spin_box<int, QSpinBox>(n_params);
		_m_spin_box = create_spin_box<int, QSpinBox>(m_params);

		_view_reset_button = new QPushButton("Reset View");
		_log_button = new QPushButton("Log Data");

		_function_checkboxes = {
			{"e^x", new QCheckBox("e^x")},
			{"sin(x)", new QCheckBox("sin(x)")},
			{"cos(x)", new QCheckBox("cos(x)")},
			{"ln(1 + x)", new QCheckBox("ln(1 + x)")},
			{"arctan(x)", new QCheckBox("arctan(x)")}
		};

		connect_controls();

		QVBoxLayout* function_layout = new QVBoxLayout();
		for (auto& [name, checkbox] : _function_checkboxes) {
			function_layout->addWidget(checkbox);
		}

		control_layout->addWidget(new QLabel("Order:"));
		control_layout->addWidget(_order_spin_box);
		control_layout->addWidget(new QLabel("Padé n:"));
		control_layout->addWidget(_n_spin_box);
		control_layout->addWidget(new QLabel("Padé m:"));
		control_layout->addWidget(_m_spin_box);
		control_layout->addLayout(function_layout);
		control_layout->addWidget(_view_reset_button);
		control_layout->addWidget(_log_button);

		_control_panel->setMinimumWidth(250);

		QHBoxLayout* main_layout = new QHBoxLayout();
		main_layout->addWidget(_plot, 1);
		main_layout->addWidget(_control_panel);

		setLayout(main_layout);
		update_plot();
	}

private:
	template <typename Number, typename SpinBoxType>
	static SpinBoxType* create_spin_box(const VariableParams<Number>& params) {
		SpinBoxType* spin_box = new SpinBoxType();
		spin_box->setRange(params.min, params.max);
		spin_box->setValue(params.initial);
		spin_box->setSingleStep(params.step);
		return spin_box;
	}

	void connect_controls() {
		connect(_order_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);
		connect(_n_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);
		connect(_m_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);
		connect(_view_reset_button, &QPushButton::clicked, this, [this]() {
			_default_axis_ranges = true;
			update_plot();
		});
		connect(_log_button, &QPushButton::clicked, this, &PlotWindow::log_data);
		for (auto& [name, checkbox] : _function_checkboxes) {
			connect(checkbox, &QCheckBox::stateChanged, this, &PlotWindow::update_plot);
		}
	}

	void update_plot() {
		const int order = _order_spin_box->value();
		const int n = _n_spin_box->value();
		const int m = _m_spin_box->value();

		_plot->clearGraphs();

		std::vector<std::pair<QString, std::vector<double>>> selected_functions;
		if (_function_checkboxes["e^x"]->isChecked())
			selected_functions.emplace_back("e^x", maclaurin_exp(order));
		if (_function_checkboxes["sin(x)"]->isChecked())
			selected_functions.emplace_back("sin(x)", maclaurin_sin(order));
		if (_function_checkboxes["cos(x)"]->isChecked())
			selected_functions.emplace_back("cos(x)", maclaurin_cos(order));
		if (_function_checkboxes["ln(1 + x)"]->isChecked())
			selected_functions.emplace_back("ln(1 + x)", maclaurin_ln(order));
		if (_function_checkboxes["arctan(x)"]->isChecked())
			selected_functions.emplace_back("arctan(x)", maclaurin_arctan(order));

		for (size_t i = 0; i < selected_functions.size(); ++i) {
			const auto& [name, coeffs] = selected_functions[i];

			QVector<double> x, y;
			for (double xi = -1.0; xi <= 1.0; xi += 0.01) {
				double yi = evaluate_polynomial(coeffs, xi);
				x.append(xi);
				y.append(yi);
			}

			_plot->addGraph();
			_plot->graph(i)->setName(name);
			_plot->graph(i)->setData(x, y);
			_plot->graph(i)->setPen(QPen(QColor::fromHsv(60 * i, 255, 255)));

			// Добавление аппроксимации Паде
			auto [num, den] = alfi::ratf::pade(coeffs, n, m);
			QVector<double> pade_y;
			for (double xi : x) {
				double yi = evaluate_rational(num, den, xi);
				pade_y.append(yi);
			}

			_plot->addGraph();
			_plot->graph(i + selected_functions.size())->setPen(QPen(Qt::DashLine));
			_plot->graph(i + selected_functions.size())->setData(x, pade_y);
		}

		_plot->replot();
	}

	void log_data() {
		std::cout << "=== Log Data ===\n";
		for (const auto& [name, checkbox] : _function_checkboxes) {
			if (checkbox->isChecked()) {
				std::cout << "Function: " << name.toStdString() << "\n";
			}
		}
	}

	static double evaluate_polynomial(const std::vector<double>& coeffs, double x) {
		double result = 0.0;
		for (double coeff : coeffs) {
			result = result * x + coeff;
		}
		return result;
	}

	static double evaluate_rational(const std::vector<double>& num, const std::vector<double>& den, double x) {
		return evaluate_polynomial(num, x) / evaluate_polynomial(den, x);
	}

private:
	bool _default_axis_ranges = true;
	QCustomPlot* _plot;
	QWidget* _control_panel;
	QSpinBox *_order_spin_box, *_n_spin_box, *_m_spin_box;
	QPushButton *_view_reset_button, *_log_button;
	std::map<QString, QCheckBox*> _function_checkboxes;
};

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	PlotWindow window;
	window.setWindowTitle("Function Approximation");
	window.resize(1200, 600);
	window.show();
	return QApplication::exec();
}