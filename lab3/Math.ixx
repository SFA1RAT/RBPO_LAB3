module;
#include <math.h>
#include <ostream>
export module Math;

export class Complex
{
private:
	double m_re;
	double m_im;
public:
	//  онструкторы дл€ инициализации комплексных чисел
	Complex(double real)
	{
		m_re = real;
		m_im = 0;
	}
	Complex()
	{
		m_re = 0;
		m_im = 0;
	}
	Complex(double real, double imag)
	{
		m_re = real;
		m_im = imag;
	}
	// —татические методы дл€ создани€ комплексных чисел из экспоненциальной и алгебраической форм
	static Complex FromExponentialForm(double _mod, double _arg)
	{
		Complex exp_obj;
		exp_obj.m_re = _mod * cos(_arg);
		exp_obj.m_im = _mod * sin(_arg);
		return exp_obj;
	}
	static Complex FromAlgebraicForm(double real, double imag)
	{
		Complex alg_obj(real, imag);
		//alg_obj.m_re = real;
		//alg_obj.m_im = imag;
		return alg_obj;
	}
	// ћетоды доступа к действительной и мнимой част€м, модулю и аргументу комплексного числа
	double Re() const
	{
		return m_re;
	}
	double Im() const
	{
		return m_im;
	}
	double Mod() const
	{
		return sqrt(m_re * m_re + m_im * m_im);
	}
	double Arg() const
	{
		return atan2(m_im, m_re);
	}
	// ќператор приведени€ типа к double
	explicit operator double() const {
		return m_re;
	}
	// ќператоры изменени€ знака, инкремента и декремента
	Complex operator-()
	{
		Complex obj(*this);
		obj.m_im *= -1;
		obj.m_re *= -1;
		return obj;
	}
	Complex& operator++()
	{
		m_re++;
		return (*this);
	}
	Complex operator++(int incr)
	{
		Complex obj(*this);
		++(*this);
		return obj;
	}
	Complex& operator--()
	{
		m_re--;
		return (*this);
	}
	Complex operator--(int decr)
	{
		Complex obj(*this);
		--(*this);
		return obj;
	}
	// ќператоры сложени€, вычитани€, умножени€ и делени€ комплексных чисел
	Complex& operator+=(Complex temp_obj)
	{
		m_re += temp_obj.m_re;
		m_im += temp_obj.m_im;
		return (*this);
	}
	Complex& operator-=(Complex temp_obj) {
		m_re -= temp_obj.m_re;
		m_im -= temp_obj.m_im;
		return (*this);
	}
	
	Complex& operator*=(Complex temp_obj) {
		double temp_re = m_re;
		double temp_im = m_im;
		m_re = temp_re * temp_obj.m_re - temp_im * temp_obj.m_im;
		m_im = temp_re * temp_obj.m_im + temp_im * temp_obj.m_re;
		return (*this);
	}
	Complex& operator/=(Complex temp_obj) {
		double temp_re1 = m_re, temp_im1 = m_im;
		double temp_re2 = temp_obj.m_re, temp_im2 = temp_obj.m_im;
		m_re = (temp_re1 * temp_re2 + temp_im1 * temp_im2) / (pow(temp_re2, 2) + pow(temp_im2, 2));
		m_im = (temp_re2 * temp_im1 - temp_re1 * temp_im2) / (pow(temp_re2, 2) + pow(temp_im2, 2));
		return (*this);
	}
	// ƒружественные операторы сложени€, вычитани€, умножени€ и делени€ комплексных чисел
	friend Complex operator+ (const Complex& temp_obj1, const Complex& temp_obj2);
	friend Complex operator- (const Complex& temp_obj1, const Complex& temp_obj2);
	friend Complex operator* (const Complex& temp_obj1, const Complex& temp_obj2);
	friend Complex operator/ (const Complex& temp_obj1, const Complex& temp_obj2);

	// ƒружественные операторы дл€ литералов i
	friend Complex operator ""i(long double imag);
	friend Complex operator ""i(unsigned long long imag);

	// ƒружественный оператор вывода в поток
	friend std::ostream& operator<<(std::ostream& stream, const Complex& temp_obj);
};
export Complex operator+(const Complex& temp_obj1, const Complex& temp_obj2)
{
	return Complex(temp_obj1.m_re + temp_obj2.m_re, temp_obj1.m_im + temp_obj2.m_im);
}
export Complex operator-(const Complex& temp_obj1, const Complex& temp_obj2)
{
	return Complex(temp_obj1.m_re - temp_obj2.m_re, temp_obj1.m_im - temp_obj2.m_im);
}
export Complex operator*(const Complex& temp_obj1, const Complex& temp_obj2)
{
	return Complex((temp_obj1.m_re * temp_obj2.m_re - temp_obj1.m_im * temp_obj2.m_im),
		(temp_obj1.m_re * temp_obj2.m_im + temp_obj1.m_im * temp_obj2.m_re));
}
export Complex operator/(const Complex& temp_obj1, const Complex& temp_obj2)
{
	return Complex((temp_obj1.m_re * temp_obj2.m_re + temp_obj1.m_im * temp_obj2.m_im) /
		(temp_obj2.m_re * temp_obj2.m_re + temp_obj2.m_im * temp_obj2.m_im),
		(temp_obj2.m_re * temp_obj1.m_im - temp_obj1.m_re * temp_obj2.m_im) /
		(temp_obj2.m_re * temp_obj2.m_re + temp_obj2.m_im * temp_obj2.m_im));
}
export Complex operator""i(long double imag)
{
	return Complex(0.0, static_cast<double>(imag));
}
export Complex operator""i(unsigned long long imag)
{
	return Complex(0.0, static_cast<double>(imag));
}
export std::ostream& operator<<(std::ostream& stream, const Complex& temp_obj)
{
	if (temp_obj.m_im < 0)
	{
		stream << temp_obj.m_re << " " << temp_obj.m_im << "i";
	}
	else
	{
		stream << temp_obj.m_re << " + " << temp_obj.m_im << "i";
	}
	return stream;
}

export int FindGreatestCommonDivisor(int a, int b)
{
	int r;
	if (a < 0)
		a *= -1;
	if (b < 0)
		b *= -1;
	while (true)
	{
		if (b == 0)
			return a;
		r = a % b;
		a = b;
		b = r;
	}
}
export int FindLeastCommonMultiple(int x, int y) {
	return abs(x * y) / FindGreatestCommonDivisor(x, y);
}

export class Rational {
	int m_nominator;
	int m_denominator;

public:
	// ћетод дл€ нормализации рационального числа
	void normalize()
	{
		int nod = FindGreatestCommonDivisor(m_nominator, m_denominator);
		m_nominator /= nod;
		m_denominator /= nod;
		if (m_denominator < 0) {
			m_denominator *= -1;
			m_nominator *= -1;
		}
	}
	//  онструкторы дл€ инициализации рациональных чисел
	Rational()
	{
		m_nominator = 0;
		m_denominator = 1;
	}
	Rational(int nom, int denom) {
		m_denominator = denom;
		m_nominator = nom;
		normalize();
	}
	Rational(int nom) {
		m_nominator = nom;
		m_denominator = 1;
	}
	// ћетоды доступа к числителю и знаменателю
	int Nominator() const {
		return m_nominator;
	}
	int Denominator() const {
		return m_denominator;
	}
	// ќператор приведени€ типа к double
	explicit operator double() const {
		return double(m_nominator) / m_denominator;
	}
	// ќператор изменени€ знака, инкремента и декремента
	Rational operator-() {
		Rational obj(*this);
		obj.m_nominator *= -1;
		return obj;
	}
	Rational& operator++ () {
		m_nominator += m_denominator;
		return (*this);
	}
	Rational operator++ (int param) {
		Rational obj(*this);
		(*this).m_nominator += m_denominator;
		return obj;
	}
	Rational& operator-- () {
		m_nominator -= m_denominator;
		return (*this);
	}
	Rational operator-- (int param) {
		Rational obj(*this);
		(*this).m_nominator -= m_denominator;
		return obj;
	}
	// ќператоры сложени€, вычитани€, умножени€ и делени€ рациональных чисел
	Rational& operator+=(Rational temp_obj) {
		int new_den = FindLeastCommonMultiple(m_denominator, temp_obj.m_denominator);
		m_nominator = new_den / m_denominator * m_nominator;
		m_nominator += new_den / temp_obj.m_denominator * temp_obj.m_nominator;
		m_denominator = new_den;
		normalize();
		return (*this);
	}
	Rational& operator-=(Rational temp_obj) {
		int new_d = FindGreatestCommonDivisor(m_denominator, temp_obj.m_denominator);
		m_nominator = new_d / m_denominator * m_nominator;
		m_nominator -= new_d / temp_obj.m_denominator * temp_obj.m_nominator;
		m_denominator = new_d;
		normalize();
		return (*this);
	}
	Rational& operator*=(Rational temp_obj) {
		m_denominator *= temp_obj.m_denominator;
		m_nominator *= temp_obj.m_nominator;
		normalize();
		return (*this);
	}
	Rational& operator/=(Rational temp_obj) {
		m_denominator *= temp_obj.m_nominator;
		m_nominator *= temp_obj.m_denominator;
		normalize();
		return (*this);
	}
	// ƒружественные операторы сложени€, вычитани€, умножени€ и делени€ рациональных чисел
	friend Rational operator+ (const Rational& temp_obj1, const Rational& temp_obj2);
	friend Rational operator- (const Rational& temp_obj1, const Rational& temp_obj2);
	friend Rational operator* (const Rational& temp_obj1, const Rational& temp_obj2);
	friend Rational operator/(const Rational& temp_obj1, const Rational& temp_obj2);

	friend bool operator==(const Rational& temp_obj1, const Rational& temp_obj2);
	friend bool operator>(const Rational& temp_obj1, const Rational& temp_obj2);
	friend bool operator<(const Rational& temp_obj1, const Rational& temp_obj2);
	friend bool operator>=(const Rational& temp_obj1, const Rational& temp_obj2);
	friend bool operator<=(const Rational& temp_obj1, const Rational& temp_obj2);

	friend std::ostream& operator<<(std::ostream& stream, const Rational& temp_obj);
};

export Rational operator+ (const Rational& temp_obj1, const Rational& temp_obj2) {
	int denominator = FindLeastCommonMultiple(temp_obj1.m_denominator, temp_obj2.m_denominator);
	int nominator = denominator / temp_obj1.m_denominator * temp_obj1.m_nominator;
	nominator += denominator / temp_obj2.m_denominator * temp_obj2.m_nominator;
	return Rational{ nominator, denominator };
}

export Rational operator-(const Rational& temp_obj1, const Rational& temp_obj2)
{
	int denominator = FindLeastCommonMultiple(temp_obj1.m_denominator, temp_obj2.m_denominator);
	int nominator = denominator / temp_obj1.m_denominator * temp_obj1.m_nominator;
	nominator -= denominator / temp_obj2.m_denominator * temp_obj2.m_nominator;
	return Rational{ nominator, denominator };
}

export Rational operator*(const Rational& temp_obj1, const Rational& temp_obj2)
{
	return Rational{ temp_obj1.m_nominator * temp_obj2.m_nominator, temp_obj2.m_denominator * temp_obj1.m_denominator };
}

export Rational operator/(const Rational& temp_obj1, const Rational& temp_obj2)
{
	return Rational{ temp_obj1.m_nominator * temp_obj2.m_denominator,temp_obj1.m_denominator * temp_obj2.m_nominator };
}

export bool operator==(const Rational& temp_obj1, const Rational& temp_obj2)
{
	return temp_obj1.m_nominator == temp_obj2.m_nominator && temp_obj1.m_denominator == temp_obj2.m_denominator;
}

export bool operator>(const Rational& temp_obj1, const Rational& temp_obj2)
{
	int den = FindLeastCommonMultiple(temp_obj1.m_denominator, temp_obj2.m_denominator);
	return den / temp_obj1.m_denominator * temp_obj1.m_nominator > den / temp_obj2.m_denominator * temp_obj2.m_nominator;
}
export bool operator<(const Rational& temp_obj1, const Rational& temp_obj2)
{
	int den = FindLeastCommonMultiple(temp_obj1.m_denominator, temp_obj2.m_denominator);
	return den / temp_obj1.m_denominator * temp_obj1.m_nominator < den / temp_obj2.m_denominator * temp_obj2.m_nominator;
}
export bool operator>=(const Rational& temp_obj1, const Rational& temp_obj2)
{
	int den = FindLeastCommonMultiple(temp_obj1.m_denominator, temp_obj2.m_denominator);
	return den / temp_obj1.m_denominator * temp_obj1.m_nominator >= den / temp_obj2.m_denominator * temp_obj2.m_nominator;
}
export bool operator<=(const Rational& temp_obj1, const Rational& temp_obj2)
{
	int den = FindLeastCommonMultiple(temp_obj1.m_denominator, temp_obj2.m_denominator);
	return den / temp_obj1.m_denominator * temp_obj1.m_nominator <= den / temp_obj2.m_denominator * temp_obj2.m_nominator;
}

export std::ostream& operator<<(std::ostream& stream, const Rational& temp_obj) {
	stream << temp_obj.m_nominator << "|" << temp_obj.m_denominator;
	return stream;
}

export Complex f(const Complex& z)
{
	Complex a(0, 0);
	Complex result = z * z * z + (1 + 2 * a) * z * z + (1 - 2 * a) * (z * z * z * z * z);
	return result;
}

export Rational f(const Rational& r)
{
	Rational a(1, 2);
	Rational temp(r * r * r * r * r);
	double chisl = temp.Nominator();
	double znam = temp.Denominator();
	double drob = chisl / znam;
	Rational result = r * r * r + (1 + 2 * a) * r * r + (1 - 2 * a) * pow(drob, 2.5);
	return result;
}

export double f(const double& d)
{
	double a = 0.5;
	double drob = (d * d * d * d * d);
	double result = d * d * d + (1 + 2 * a) * d * d + (1 - 2 * a) * pow(drob, 2.5);
	return result;
}