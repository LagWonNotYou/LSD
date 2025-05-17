//	L.S.D. (Lhamo-Sergi's Dependencies) //DA FINIRE Ridefinire funzioni, cambiare i nomi in lowercase
/*	INDEX 
	(0)Base
	(1)Source Code
	(2)Math
	(3)Vectors
	(4)Data Structs
	
//(0) BASE CONCEPTS AND EXPLENATION
	explenation of T class, vector<T>, ecc
	
//(1) SOURCE CODE	

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <utility>
#include <fstream>
	
#include "LSD.h" //v1

/*  Function List

	(0) Base
	void help
	void function_list
	
	(1) Math
	bool is_prime
	bool is_even
	bool p_sqr
	bool mult_of
	int abs_val
	int fact
	int gap
	int rand_btw
	int max_btw
	int min_btw 
	
	(2) Vectors
	void fill_vect
	void rand_fill_vect
	void print_vect
	void shift_vect
	bool check_vect
	int max_vect 
	int min_vect
	bool check_vect_pos
	int max_vect_pos 
	int min_vect_pos
	int count_equal_vect
	int count_max_vect 
	int count_min_vect 
	int sum_vect
	float AverageVect 
	void ParallelSum 
	void ParallelSub
	void IncSelSort
	void DecSelSort 
	void ParallelDecSelSort

	(3) Data Structs
	vector_plus<T>:
		vector<T> find_all
		int find_first
		int fast_find
		void sort
		bool is_sorted
		operator>> (cin)
		operator<< (cout)
	matrix<T>:
		void randFill
		void edit_height
		void edit_length
		vector<T> row
		vector<T> column
		void clear
		void reset
		bool is_square
		void set
		T get
		T determinant
		T sum
		T max
		T min
		int max_pos
		int min_pos
		int* data
		vector<T>::iterator begin
		vector<T>::iterator end
		operator= 
		operator<< (cout)
*/

//(2) MATH FUNCTIONS

bool is_prime(int Number_to_check) {
	int interval, multiples=0;
	for(interval=1;interval<=Number_to_check;interval++) {
		if(Number_to_check%interval==0) multiples++;
	}
	return (multiples == 2);
}

int abs_val(int Number) {
	return abs(Number);
}
			
bool is_even(int n) {
	return (n % 2 == 0);
}
	
int fact(int Number) {
	int interval, result=1;
	for(interval=2;interval<=Number;interval++) {
	result=result*interval;	
	}
		return result;
}
	
int rand_btw(int Random_Number_1,int Random_Number_2) {
	if(Random_Number_1>Random_Number_2) {
		std::swap(Random_Number_1,Random_Number_2);
	}
		return rand()%(Random_Number_2-Random_Number_1+1)+Random_Number_1;
}
	
int gap(int Number_1, int Number_2, int Number_3) {
	return std::max(std::max(Number_1,Number_2),Number_3)-std::min(std::min(Number_1,Number_2),Number_3);
}

bool p_sqr(int n) {
	int m = sqrt(n);
	return (m * m == n);
}

bool mult_of(int n, int mult) {
	return (n % mult == 0);
}

//(3) VECTOR FUNCTIONS

void FillVect(int Vector[], int Length) {
	int i;
	for(i=0;i<Length;i++) {
		std::cin>>Vector[i];
	}
}

void RandFillVect(int Vector[], int Length, int Random_Number_1, int Random_Number_2) {
	int interval;
	if(Random_Number_1>Random_Number_2) {
		std::swap(Random_Number_1,Random_Number_2);
	}
	for(interval=0;interval<Length;interval++) {
		Vector[interval]=rand_btw(Random_Number_1,Random_Number_2);
	}
}

void PrintVect(int Vector[], int Length) {
	int interval;
	std::cout<<std::endl;
	for(interval=0;interval<Length;interval++) {
		std::cout<<"\t"<<Vector[interval];
	}
	std::cout<<std::endl;
}

void ShiftVect(int Vector[], int Length, int Shift_Amount) {
	int temp_vector[Length], interval;
	if(Shift_Amount>0) {
		for(interval=0;interval<Shift_Amount;interval++) {
			temp_vector[interval]=Vector[interval-Shift_Amount+Length];
		}
		for(interval=Shift_Amount;interval<Length;interval++) {
			temp_vector[interval]=Vector[interval-Shift_Amount];
		}
	}
	else {
		for(interval=0;interval<Length;interval++) {
			temp_vector[interval]=Vector[(interval-Shift_Amount)%Length];
		}
	}
		for(interval=0;interval<Length;interval++) {
			Vector[interval]=temp_vector[interval];
		}
}
	
bool CheckVect(int Vector[], int Length, int Number_to_check) {
	int interval;
	for(interval=0;interval<Length;interval++)
	{
		if(Vector[interval]==Number_to_check)
			return true;
	}
		return false;
}

int MaxVect(int Vector[], int Length) {
	int interval, max;
	max=Vector[0];
	for(interval=1;interval<Length;interval++) {
		if(Vector[interval]>max)
		max=Vector[interval];
	}
		return max;
}

int MinVect(int Vector[], int Length) {
	int interval, min;
		min=Vector[0];
		for(interval=1;interval<Length;interval++) {
			if(Vector[interval]<min)
			min=Vector[interval];
		}
			return min;
}
	
int CheckVectPos(int Vector[], int Length, int Number_to_check) {
	int interval, position=0;
	for(interval=0;interval<Length;interval++) {
		if(Vector[interval]==Number_to_check)
			position=interval;
	}
			return position;
}

int MaxVectPos(int Vector[], int Length) {
	int interval, max, max_position=0;
	max=Vector[0];
	for(interval=1;interval<Length;interval++) {
		if(Vector[interval]>max) {
			max=Vector[interval];
			max_position=interval;
		}
	}
		return max_position;
}

int MinVectPos(int Vector[], int Length) {
	int interval, min, min_position=0;
	min=Vector[0];
	for(interval=1;interval<Length;interval++) {
		if(Vector[interval]<min) {
			min=Vector[interval];
			min_position=interval;
		}
	}
		return min_position;
}

int CountEqualVect(int Vector[], int Length, int Number_to_check) {
	int interval, counter=0;
	for(interval=0;interval<Length;interval++) {
		if(Vector[interval]==Number_to_check)
			counter++;
	}
		return counter;
}

int CountMaxVect(int Vector[], int Length, int Number_to_check) {
	int interval, counter=0;
	for(interval=0;interval<Length;interval++) {
		if(Vector[interval]>Number_to_check)
			counter++;
	}
		return counter;
}

int CountMinVect(int Vector[], int Length, int Number_to_check) {
	int interval, counter=0;
	for(interval=0;interval<Length;interval++) {
		if(Vector[interval]<Number_to_check)
			counter++;
	}
	return counter;
}

int SumVect(int Vector[], int Length) {
	int counter=0, interval;
	for(interval=0;interval<Length;interval++) {
		counter=counter+Vector[interval];
	}
		return counter;
}
	
float AverageVect(int Vector[], int Length) {
	int interval; 
	float counter=0, result;
	for(interval=0;interval<Length;interval++) {
		counter=counter+Vector[interval];
	}
	result=counter/Length;
		return result;
}

void ParallelSum(int First_Vector[], int Second_Vector[], int Result_Vector[], int Length) {
	int interval;
	for(interval=0;interval<Length;interval++) {
		Result_Vector[interval]=First_Vector[interval]+Second_Vector[interval];
	}
}

void ParallelSub(int First_Vector[], int Second_Vector[], int Result_Vector[], int Length) {
	int interval;
	for(interval=0;interval<Length;interval++) {
		Result_Vector[interval]=First_Vector[interval]-Second_Vector[interval];
	}
}

void IncSelSort(int Vector[], int Length) {
	int i, j, min, pos_min;
	for(i=0;i<Length-1;i++) {														
		min=Vector[i];
		pos_min=i;
		for(j=i+1;j<Length;j++) {											
			if(Vector[j]<min){
				min=Vector[j];
				pos_min=j;
			}
		}
		if(pos_min!=i)									
			std::swap(Vector[pos_min],Vector[i]);
	}
}

void DecSelSort(int Vector[], int Length) {
	int i, j, max, pos_max;
	for(i=0;i<Length-1;i++) {					
		max=Vector[i];
		pos_max=i;
		for(j=i+1;j<Length;j++) {				
			if(Vector[j]>max){
				max=Vector[j];
				pos_max=j;
			}
		}
		if(pos_max!=i)
			std::swap(Vector[pos_max],Vector[i]);
	}
}

void ParallelDecSelSort(int Independent_Vector[], int Dependent_Vector[], int Length) {
	int i, j, max, pos_max;
	for(i=0;i<Length-1;i++) {					
		max=Independent_Vector[i];
		pos_max=i;
		for(j=i+1;j<Length;j++) {				
			if(Independent_Vector[j]>max) {
				max=Independent_Vector[j];
				pos_max=j;
			}
		}
		if(pos_max!=i) {
			std::swap(Independent_Vector[pos_max],Independent_Vector[i]);
			std::swap(Dependent_Vector[pos_max],Dependent_Vector[i]);				
		}
	}
}

//(4) DATA STRUCTS

#ifdef _GLIBCXX_VECTOR
#ifndef DATA_STRUCTS_GUARD
#define DATA_STRUCTS_GUARD

template<typename T> using func = bool (T, T);

template<typename T> bool default_comp(T lhs, T rhs) {return lhs < rhs;}

template<class T> class vector_plus : public std::vector<T> {
	template<class C> friend std::ostream& operator<<(std::ostream&, const vector_plus<C>&);
	template<class C> friend std::istream& operator>>(std::istream&, vector_plus<C>&);
	public:
		vector_plus<T>() {
		}
		vector_plus<T>(std::initializer_list<T> init) {
			for (typename std::initializer_list<T>::iterator i = init.begin(); i < init.end(); i++) {
				this->push_back(*i);
			}
		}		
		vector_plus<T>(std::vector<T> vc) {
			for (typename std::vector<T>::iterator i = vc.begin(); i < vc.end(); i++) {
				this->push_back(*i);
			}
		}
		int find_first(T value) {
			for (int i = 0; i < this->size(); i++) {
				if (this->at(i) == value) return i;
			}
			return -1;
		}
		int fast_find(T value) {
			int l = 0, r = this->size() - 1, m;
			while (l <= r) {
				m = std::floor((l + r)/2);
				if (this->at(m) > value) r = m - 1;
				else if (this->at(m) < value) l = m + 1;
				else return m;
			}
			return -1;
		}
		std::vector<T> find_all(T value) {
			std::vector<T> out;
			for (int i = 0; i < this->size(); i++) {
				if (this->at(i) == value) out.push_back(i);
			}
			return out;
		}
		void sort(func<T> comp = default_comp<T>) {
			#ifndef _GLIBCXX_ALGORITHM
			throw std::runtime_error("Missing library: <algorithm>");
			#endif
			#ifdef _GLIBCXX_ALGORITHM
			std::sort(this->begin(), this->end(), comp);
			#endif
		}
		bool is_sorted(func<T> comp = default_comp<T>) {
			for (int i = 0; i < this->size() - 1; i++) {
				if(!(comp(this->at(i), this->at(i + 1)))) return false;
			}
			return true;
		}
};

template<typename T> class matrix {
	private:
		int height, length;
		std::vector<T> storage;
	public:
		matrix(int h, int l) {
			length = l;
			height = h;
			for (int i = 0; i < h * l; i++) storage.push_back(0);
		}
		matrix(int h, int l, std::initializer_list<T> init) {
			length = l;
			height = h;
			if (init.size() == (h * l)) {
				for (T t : init) {
					storage.push_back(t);
				}
			}
			else throw std::invalid_argument("Matrix initialization failed: size error");
		}
		matrix& operator=(const std::initializer_list<T>& init) {
			if (init.size() == storage.size()) {
				int i = 0;
				for (typename std::initializer_list<T>::iterator it = init.begin(); it != init.end(); it++) {
					storage.at(i) == *it;
					i++;
				}
			}
			else {
				throw std::invalid_argument("Matrix assignment failed: size error");
			}
			return *this;
		}
		void randFill(int lb = -1, int ub = 1) {
			#ifdef _GLIBCXX_CSTDLIB
			for (typename std::vector<T>::iterator it = storage.begin(); it != storage.end(); it++) {
				*it = (rand() % (ub - lb + 1)) + lb;
			}
			#endif
			#ifndef _GLIBCXX_CSTDLIB
			throw std::runtime_error("Missing library: <cstdlib>");
			#endif
		}
		template<typename C> friend std::ostream& operator<< (std::ostream&, matrix<C>&);
		T get(int h, int l) {return storage.at((h * length) + l);}
		void set(int h, int l, T value) {storage.at((h * length) + l) = value;}
		void edit_height(int h) {height = h;}
		void edit_length(int l) {length = l;}
		void clear() {storage.clear();}
		void reset() {for (int i = 0; i < storage.size(); i++) storage[i] = 0;}
		bool is_square() {return length == height;}
		std::vector<T> row(int req) {
			std::vector<T> out;
			for (int i = 0; i < length; i++) {
				out.push_back(storage[(req * length) + i]);
			}
			return out;
		}
		std::vector<T> column(int req) {
			std::vector<T> out;
			for (int i = 0; i < height; i++) {
				out.push_back(storage[(length * i) + req]);
			}
			return out;
		}
		T determinant() {
			if (!is_square()) throw std::runtime_error("Non-square matrix has no determinant");
			if (height > 3) throw std::runtime_error("Nuh uh (determinant is available up to 3x3)");
			switch (height) {
				case 1:
					return storage[0];
				case 2:
					return (storage[0] * storage[3]) - (storage[1] * storage[2]);
				case 3:
					return (storage[0] * storage[4] * storage[8]) + (storage[1] * storage[5] * storage[6]) + (storage[2] * storage[3] * storage[7]) - (storage[2] * storage[4] * storage[6]) - (storage[0] * storage[5] * storage[7]) - (storage[1] * storage[3] * storage[8]);
			}
			return -1;
		}
		T sum() {
			int out = 0;
			for (int t : storage) out += t;
			return out;
		}
		T max() {
			int m = 0;
			for (int t : storage) if (t > m) m = t;
			return m;
		}
		T min() {
			int m = storage[0];
			for (int t : storage) if (t < m) m = t;
			return m;
		}
		int max_pos() {
			int pos = 1, m = storage[0];
			for (int i = 1; i < storage.size(); i++) {
				if (storage[i] > m) {
					pos = i;
					m = storage[i];
				}
			}
			return pos;
		}
		int min_pos() {
			int pos = 1, m = storage[0];
			for (int i = 1; i < storage.size(); i++) {
				if (storage[i] < m) {
					pos = i;
					m = storage[i];
				}
			}
			return pos;
		}
		int* data() {
			return storage.data();
		}
		typename std::vector<T>::iterator begin() {
			return storage.begin();
		}
		typename std::vector<T>::iterator end() {
			return storage.end();
		}
};

template<class C> std::ostream& operator<<(std::ostream& os, const vector_plus<C>& vc) {
	for (int i = 0; i < vc.size(); i++) {
		os << vc.at(i) << " ";
	}
	return os;
}

template<class C> std::istream& operator>>(std::istream& is, vector_plus<C>& vc) {
	C temp;
	is >> temp;
	vc.push_back(temp);
	return is;
}

template<typename C> std::ostream& operator<< (std::ostream& os, matrix<C>& mx) {
	for (int i = 1; i <= mx.storage.size(); i++) {
		os << mx.storage[i - 1];
		if (i % mx.length == 0) os << std::endl;
		else os << ' ';
	}
	return os;
}
#endif
#endif
