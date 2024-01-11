#include <iostream>
import Math;
using namespace std;

int main() {
    setlocale(LC_ALL, "Russian");

    double re, im;
    cout << "Введите действительную и мнимую часть комплексного числа:" << endl;
    cin >> re >> im;
    Complex c_num(re, im);

    double chisl, znam;
    cout << "Введите числитель и знаменатель для рационального числа:" << endl;
    cin >> chisl >> znam;
    Rational r_num(chisl, znam);

    double d_num;
    cout << "Введите вещественное число:" << endl;
    cin >> d_num;

    Rational answ = f(r_num);
    double a = answ.Nominator(), b = answ.Denominator();

    cout << "\nРезультаты вычислений:\n";
    cout << "-------------------------\n";
    cout << "f(" << c_num << ") = " << f(c_num) << "  (комплексное число)\n";
    cout << "f(" << r_num << ") = " << f(r_num) << " = " << a << " / " << b << "  (рациональное число)\n";
    cout << "f(" << d_num << ") = " << f(d_num) << "  (вещественное число)\n";

    return 0;
}