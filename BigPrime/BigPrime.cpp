/*
	附加包含目录:gmp
	附加库目录:MinGW\lib
	附加依赖项:libgmp.dll.a;libgmpxx.dll.a
*/
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
using namespace std;

const static int maxLength = 120;
const static int primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};

class BigPrime {
	mpz_t n;						//存储大素数
	gmp_randstate_t grt;			//随机参数
	int t;							//精度参数 1-(1/4)^t
	int Length;						//二进制长度	
	void RandSet();					//将n随机化
	bool isPrime();					//判断是否为素数
public:
	void copy(mpz_t& out);
	BigPrime(int t, int Length);
	~BigPrime();
};

void BigPrime::copy(mpz_t& out) {
	mpz_add_ui(out, n, 0);
}

BigPrime::~BigPrime() {
	gmp_randclear(grt);
	mpz_clear(n);
}

bool BigPrime::isPrime() {
	int s = 0, flag = 1;
	mpz_t r1, r2, b, n_sub_one, k, r;
	mpz_init(r1);
	mpz_init(r2);
	mpz_init(b);
	mpz_init(n_sub_one);
	mpz_init(k);
	mpz_init(r);

	mpz_cdiv_r_ui(r, n, 2);
	if (mpz_cmp_ui(r, 0) == 0) return false;	//偶数直接返回false
	
	/*	2^s*k = n-1	*/
	mpz_sub_ui(k, n, 1);
	mpz_cdiv_r_ui(r, k, 2);
	while (mpz_cmp_ui(r, 0) == 0) {
		s++;
		mpz_cdiv_q_ui(k, k, 2);
		mpz_cdiv_r_ui(r, k, 2);
	}	

	for (int i = 0; i < t; i++) {
		mpz_init_set_ui(b, (unsigned int)primes[rand() % 10]);

		mpz_powm(r1, b, k, n);
		mpz_sub_ui(n_sub_one, n, 1);

		if (mpz_cmp(r1, n_sub_one) == 0 || mpz_cmp_ui(r1, 1) == 0) { continue; }
		for (; i < s - 1; i++) {
			mpz_powm_ui(r2, r1, 2, n);
			if (mpz_cmp_ui(r2, 1) == 0) flag = 0;
			mpz_add_ui(r1, r2, 0);
		}
		if (mpz_cmp_ui(r2, 1) != 0)  flag = 0;
	}
	mpz_clear(r1);
	mpz_clear(r2);
	mpz_clear(b);
	mpz_clear(n_sub_one);
	mpz_clear(k);
	mpz_clear(r);
	return flag;
}	

void BigPrime::RandSet() {
	mpz_urandomb(n, grt, Length);
}

BigPrime::BigPrime(int t, int Length) {
	this->t = t;
	this->Length = Length;
	gmp_randinit_default(grt);			//设置随机数生成算法为默认
	gmp_randseed_ui(grt, clock());		//设置随机化种子
	mpz_init(n);						//初始化n为0
	while(1){
		RandSet();						//初始化n为二进制长度为Length的随机数
		if (isPrime()) {
			gmp_printf("%Zd\n", n);
			break;
		}
	}
}

int main()
{
	double T = 0;
	int COUNT = 1000;
	for (int i = 0; i < COUNT; i++) {
		mpz_t p2;
		mpz_init(p2);
		BigPrime p = BigPrime(10, 120);
		p.copy(p2);
		if(mpz_probab_prime_p(p2, 15) > 0) T++;
	}
	printf("Accurancy: %.2f\n", 100.0 * (T / (double)COUNT));
}

