#define _USE_MATH_DEFINES
#include <random>
#include <iostream>
#include "math.h"
#include "mpfr.h"
#include <time.h>
#include <fstream>
using namespace std;


std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0.0, 1.0);
ofstream out("results.txt");

float find_qi_method1(int m)
{
    float sum = 0;
    float wi = dis(gen);
    float etta = 1 / wi - 1;
    for (int i = 0; i < m; i++)
    {
        wi = dis(gen);
        float xi = -log(wi);
        sum += xi;
    }
    if (etta > sum)
        return 1;
    else
        return 0;
}

float find_qi_method2(int m)
{
    float sum = 0;
    float wi = dis(gen);
    float etta = 1 / wi - 1;
    for (int i = 0; i < m; i++)
    {
        wi = dis(gen);
        float xi = -log(wi);
        sum += xi;
    }
    return 1 / (1 + sum);
}

float find_qi_method3(int m)
{
    mpfr_t sum;
    mpfr_init_set_d(sum, 0, MPFR_RNDN);
    float wi = dis(gen);
    //cout << endl << "              wi = " << wi << endl;
    float etta = 1 / wi - 1;
    //cout << "              etta = " << etta << endl;
    mpfr_t Etta;
    mpfr_init_set_d(Etta, etta, MPFR_RNDN);
    mpfr_t s;
    mpfr_init_set_d(s, 1, MPFR_RNDN);
    for (int i = 1; i < m; i++)
    {
        mpfr_t I;
        mpfr_init_set_d(I, i, MPFR_RNDN);
        mpfr_mul(s, s, Etta, MPFR_RNDN);
        mpfr_div(s, s, I, MPFR_RNDN);
        mpfr_clear(I);
    }
    int i = m;
    while (mpfr_get_d(s, MPFR_RNDN)>0.01)
    {
        mpfr_t I;
        mpfr_init_set_d(I, i, MPFR_RNDN);
        mpfr_mul(s, s, Etta, MPFR_RNDN);
        mpfr_div(s, s, I, MPFR_RNDN);
        mpfr_clear(I);
        mpfr_add(sum, sum, s, MPFR_RNDN);
        //cout << mpfr_get_d(s, MPFR_RNDN) << "; " << mpfr_get_d(sum, MPFR_RNDN) << "     ";
        i++;
    }
    mpfr_clear(s);
    mpfr_t exp;
    mpfr_init_set_d(exp, M_E, MPFR_RNDN);
    mpfr_pow(Etta, exp, Etta, MPFR_RNDN);
    //cout << endl << mpfr_get_d(Etta, MPFR_RNDN) << "     ";
    mpfr_div(sum, sum, Etta, MPFR_RNDN);
    //cout << mpfr_get_d(sum, MPFR_RNDN) << "     ";
    wi = mpfr_get_d(sum, MPFR_RNDN);
    //cout << 1-wi << endl;
    mpfr_clear(sum);
    mpfr_clear(exp);
    return wi;
}

float find_qi_method4(int m)
{
    float sum = 0;
    float wi;
    for (int i = 0; i < m-1; i++)
    {
        wi = dis(gen);
        float xi = -log(wi);
        sum += xi;
    }
    return sum / ((1 + sum) * (m - 1));
}

void algorithm(int m, int meth, int n)
{
    float Qin = 0;
    float Sigma = 0;
    float sig = 0;
    int stopper = 0;
    float s = 0;
    float s2 = 0;
    float qi = 0;
    int n1 = 2 * n;
    while (stopper!=1 && n<n1)
    {
        s = 0;
        s2 = 0;
        clock_t start, end;
        start = clock();
        for (int i = 0; i < n; i++)
        {
            if (meth == 1)
                qi = find_qi_method1(m);
            if (meth == 2)
                qi = find_qi_method2(m);
            if (meth == 3)
                qi = find_qi_method3(m);
            if (meth == 4)
                qi = find_qi_method4(m);
            s += qi;
            s2 += qi * qi;
        }
        end = clock();
        start = clock();
        Qin = s / n;
        Sigma = (s2 - Qin * s) / (n - 1);
        end = clock();
        cout << "Qin = " << Qin << "         Sigma = " << Sigma << endl;
        cout << "criteria:     abs(Sigma-sigma0) = " << Sigma - sig << endl;
        if (abs(Sigma - sig) <= 0.01)
        {
            float x = n + 1;
            if (Qin != 0)
                x = (66306.25 * Sigma) / (Qin * Qin);

            cout << "criteria2:     x = " << x << "      n = " << n << endl << endl;
            if (n >= x)
            {
                stopper = 1;
                cout << "Qin = " << Qin << "         Sigma = " << abs(Sigma) << endl;
                cout << "x = " << x << "      n = " << n << endl << endl;
            }
            else
                n ++;
        }
        else
        {
            sig = Sigma;
            n++;
        }
    }
    out.open("results.txt", std::ios::app);
    cout << endl << endl << "RESULT FOR      m=" << m << "   method=" << meth << endl;
    out << "\n\n RESULT FOR        m=" << m << "      method=" << meth << "\n";
    cout << "Qin = " << Qin << "         Sigma = " << abs(Sigma) << "         n = " << n << endl << endl;
    out << "Qin = " << Qin << "         Sigma = " << abs(Sigma) << "         n = " << n << "\n\n";
    out.close();

    out.open("res.csv", std::ios::app);
    out << m << ";" << meth << ";" << Qin << ";" << abs(Sigma) << ";" << n << "\n";
    out.close();
}

void algorithm2(int m, int meth, int n, int b)
{
    int iter = 0;
    float Qin = 0;
    //int b = 1000;
    float Sigma = 0;
    float sig = 0;
    int stopper = 0;
    float s = 0;
    float s2 = 0;
    float qi = 0;
    int n1 = 2 * n;
    while (stopper != 1 && n < n1)
    {
        iter++;
        s = 0;
        s2 = 0;
        clock_t start, end;
        start = clock();
        for (int i = 0; i < n; i++)
        {
            if (meth == 1)
                qi = find_qi_method1(m);
            if (meth == 2)
                qi = find_qi_method2(m);
            if (meth == 3)
                qi = find_qi_method3(m);
            if (meth == 4)
                qi = find_qi_method4(m);
            s += qi;
            s2 += qi * qi;
        }
        end = clock();
        /*if (n%100==0)
        {
            cout << "time qi x n generation   " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl;
        }*/
        start = clock();
        Qin = s / n;
        Sigma = s2 / (n - 1) + (Qin * Qin) * (n / (n - 1));
        end = clock();
        /*if (n % 100 == 0)
        {
            cout << "Qin+Sigma  time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl;
        }*/
        //cout << "Qin = " << Qin << "     Sigma = " << Sigma << "     s = " << s << endl;
        //if (n % 1000 == 0)
        //{
        cout << "Qin = " << Qin << "         Sigma = " << Sigma << endl;
        cout << "criteria:     abs(Sigma-sigma0) = " << Sigma - sig << endl;
        //}
        if (abs(Sigma - sig) <= 0.01)
        {
            float x = n + 1;
            if (Qin != 0)
                x = (66306.25 * Sigma) / (Qin * Qin);

            //if (n % 1000 == 0)
            //{
                cout << "criteria2:     x = " << x << "      n = " << n << endl << endl;
            //}
            if (n >= x)
            {
                if (b>=10)
                {
                    n -= b;
                    b /= 10;
                }
                else
                    stopper = 1;
                cout << "b = " << b << "          Qin = " << Qin << "         Sigma = " << Sigma << endl;
                cout << "x = " << x << "      n = " << n << endl << endl;
            }
            else
                n+=b;
        }
        else
        {
            sig = Sigma;
            n+=b;
        }
    }
    out.open("results.txt", std::ios::app);
    out << "\n\n RESULT FOR        m=" << m << "      method=" << meth << "\n";
    out << "Qin = " << Qin << "         Sigma = " << Sigma << "         n = " << n << "\n";
    out << "number of iterations = " << iter << "\n\n";
    out.close();
    cout << endl << endl << "RESULT FOR      m=" << m << "   method=" << meth << endl;
    cout << "Qin = " << Qin << "         Sigma = " << Sigma << "         n = " << n << endl << endl;
    cout << "number of iterations = " << iter << endl << endl;
}

void algorithm3(int m, int meth, int n, int b)
{
    int iter = 0;
    float Qin = 0;
    int Bstart = b;
    float Sigma = 0;
    float sig = 0;
    int stopper = 0;
    float s = 0;
    float s2 = 0;
    float S = 0;
    float S2 = 0;
    float qi = 0;
    int N = 0;
    int Nstart = n;
    int n1 = 2 * n;
    while (stopper != 1 && n < n1)
    {
        iter++;
        S = 0;
        S2 = 0;
        clock_t start, end;
        start = clock();
        for (int i = 0; i < n - N; i++)
        {
            if (meth == 1)
                qi = find_qi_method1(m);
            if (meth == 2)
                qi = find_qi_method2(m);
            if (meth == 3)
                qi = find_qi_method3(m);
            if (meth == 4)
                qi = find_qi_method4(m);
            s += qi;
            s2 += qi * qi;
            S += qi;
            S2 += qi * qi;
        }
        end = clock();
        start = clock();
        Qin = s / n;
        Sigma = (s2 - Qin * s) / (n - 1);
        end = clock();
        cout << "Qin = " << Qin << "         Sigma = " << Sigma << endl;
        cout << "criteria:     abs(Sigma-sigma0) = " << Sigma - sig << endl;
        if (abs(Sigma - sig) <= 0.01)
        {
            float x = n + 1;
            if (Qin != 0)
                x = (66306.25 * Sigma) / (Qin * Qin);

            cout << "criteria2:     x = " << x << "      n = " << n << endl << endl;
            if (n >= x)
            {
                if (b >= 10)
                {
                    n -= b;
                    b /= 10;
                    N = n;
                    n += b;
                    s -= S;
                    s2 -= S2;
                }
                else
                    stopper = 1;
                cout << "b = " << b << "          Qin = " << Qin << "         Sigma = " << Sigma << endl;
                cout << "x = " << x << "      n = " << n << endl << endl;
            }
            else
            {
                N = n;
                n += b;
            }
        }
        else
        {
            N = n;
            n += b;
        }
        sig = Sigma;
    }
    cout << endl << endl << "RESULT FOR      m=" << m << "   method=" << meth << endl;
    cout << "Qin = " << Qin << "         Sigma = " << Sigma << "         n = " << n << endl;
    cout << "number of iterations = " << iter << endl << endl;
    out.open("results.txt", std::ios::app);
    out << "\n\n RESULT FOR        m=" << m << "      method=" << meth << "\n";
    out << "Qin = " << Qin << "         Sigma = " << Sigma << "         n = " << n << "\n";
    out << "number of iterations = " << iter << "\n";
    out.close();

    out.open("res.csv", std::ios::app);
    out << m << ";" << meth << ";" << Qin << ";" << Sigma << ";" << n << "\n";
    out.close();
}


int main()
{
    clock_t start, end;
    out.open("results.txt");
    out << "Results and time\n\n";
    out.close();

    out.open("res.csv");
    out << "m;meth;Q;S;n\n";
    out.close();

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(1, 1, 44000, 1000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=1, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=1, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(1, 1, 170000);
    cout << "-------------------------------------------" << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(1, 2, 8000, 100);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=1, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=1, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(1, 2, 140000);
    cout << "-------------------------------------------" << endl << endl;

    start = clock();
    cout << "algorithm3" << endl;
    algorithm3(1, 3, 20000, 100);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=1, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=1, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(1, 3, 150000);

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(10, 1, 590000, 10000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=10, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(10, 2, 5000, 100);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=10, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(10, 3, 490000, 10000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=10, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm(10, 4, 2);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=10, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    cout << "===========================================" << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(100, 1, 6400000, 10000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=100, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=100, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    //algorithm(100, 1, 7000000);
    cout << "-------------------------------------------" << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm3(100, 2, 600, 10);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=100, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "=================================================================================\n\n";
    out.close();
    cout << "alg3, m=100, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(100, 2, 440000);
    cout << "-------------------------------------------" << endl << endl;

    start = clock();
    cout << "algorithm3" << endl;
    algorithm3(100, 3, 6000000, 100000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=100, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=100, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(100, 3, 1);  
    cout << "-------------------------------------------" << endl << endl;

    cout << "algorithm3" << endl;
    start = clock();
    algorithm(100, 4, 2);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=100, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=100, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(100, 4, 130000);

    cout << "===========================================" << endl << endl;
    //algorithm(10000, 1, 660000000);
    cout << "-------------------------------------------" << endl << endl;

    cout << "algorithm3" << endl;
    algorithm3(1000, 2, 60, 1);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=1000, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=1000, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    cout << "algorithm3" << endl;
    algorithm3(1000, 4, 2, 1);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=1000, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=1000, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    cout << "algorithm3" << endl;
    algorithm3(10000, 2, 2, 1);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10000, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=10000, meth=2     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    //algorithm(10000, 2, 9100000);
    cout << "-------------------------------------------" << endl << endl;

    start = clock();
    /*cout << "algorithm3" << endl;
    algorithm3(10000, 3, 6000000, 100000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10000, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out << "*********************************************************************************\n\n\n";
    out.close();
    cout << "alg3, m=10000, meth=3     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;*/
    //algorithm(10000, 3, 1);     ???????????????????????????????????????????????????????????
    cout << "-------------------------------------------" << endl << endl;
    
    cout << "algorithm3" << endl;
    start = clock();
    algorithm(10000, 4, 2);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=10000, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out.close();
    cout << "alg3, m=10000, meth=4     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;

    /*cout << "algorithm3" << endl;
    start = clock();
    algorithm3(1000, 1, 66000000, 100000);
    end = clock();
    out.open("results.txt", std::ios::app);
    out << "alg3, m=1000, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << "\n\n";
    out.close();
    cout << "alg3, m=1000, meth=1     time    " << ((double)end - start) / ((double)CLOCKS_PER_SEC) << endl << endl;
    //algorithm(10000, 4, 130000);
    cout << endl;*/
}