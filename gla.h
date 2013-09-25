#ifndef GLA_H
#define GLA_H
#include <sys/types.h>
#include <errno.h>

#define PrintErr(x , y) \
	do {\
		printf("%s\n" , (x)) ;\
		errno = (y) ;\
		exit(-1) ;\
	} while(0)

#define Print(x) printf("%s\n" , (x)) 

#define GLA_EINVAL	-1
#define GLA_ENOMEM	-2

typedef double type ;

typedef struct _gla_determinant {
	size_t n ;//行列式的阶数
	type * dat ;//行列式的元素数据
} gla_determinant ;

gla_determinant * gla_determinant_alloc(const size_t n) ;
void gla_determinant_set(gla_determinant * d , const size_t i , const size_t j , const type elem) ;
gla_determinant * gla_determinant_cpy(gla_determinant * c , const gla_determinant * d) ;
type gla_determinant_get(const gla_determinant * d , const size_t i , const size_t j) ;
void gla_determinant_free(gla_determinant * d) ;
gla_determinant * gla_determinant_cpy(gla_determinant * c , const gla_determinant * d) ;
void gla_determinant_multiplication(gla_determinant * d , const type k , const size_t i) ;
void gla_determinant_addition(gla_determinant * d , const size_t des , const size_t src) ;
void gla_determinant_swap(gla_determinant * d , const size_t i , const size_t i2) ;

type gla_determinant_alge_complement(const gla_determinant * d , const size_t i , const size_t j) ;
type gla_determinant_value(const gla_determinant * d) ;
int gla_determinant_is_diagonal(const gla_determinant * d) ;
int gla_determinant_is_upper_triangular(const gla_determinant * d) ;
int gla_determinant_is_lower_triangular(const gla_determinant * d) ;
void gla_determinant_transpose(const gla_determinant * d) ;

typedef struct _gla_vector {
	size_t n ;
	type * dat ;
} gla_vector ;

gla_vector * gla_vector_alloc(const size_t n) ;
void gla_vector_set(gla_vector * v , const size_t j , type value) ;
type gla_vector_get(const gla_vector * v , const size_t j) ;
void gla_vector_free(gla_vector * v) ;
gla_vector * gla_vector_cpy(gla_vector * des , const gla_vector * src) ;
void gla_vector_multiplication(gla_vector * v , const type k) ;
void gla_vector_addition(gla_vector * des , const gla_vector * src) ;
void gla_vector_swap(gla_vector * v1 , gla_vector * v2) ;
type gla_vector_dot_multiplication(const gla_vector * v1 , const gla_vector * v2) ;
gla_vector * gla_vector_cross_multiplication(const gla_vector * v1 , const gla_vector * v2 , gla_vector * des) ;

typedef struct _gla_matrix {
	size_t n , m ;//矩阵的行列数目
	type * dat ;//矩阵的元素数据
} gla_matrix ;

gla_matrix * gla_matrix_alloc(const size_t n , const size_t m) ;
void gla_matrix_set(gla_matrix * m , const size_t i , const size_t j , const type ele) ;
type gla_matrix_get(const gla_matrix * m , const size_t i , const size_t j) ;
void gla_matrix_free(gla_matrix * m) ;
gla_matrix * gla_matrix_cpy(gla_matrix * des , const gla_matrix * src) ;
void gla_matrix_addition(gla_matrix * des , const gla_matrix * src) ;
void gla_matrix_multiplication(gla_matrix * m , const type k) ;
gla_matrix * gla_matrix_matrix_multiplication(gla_matrix *des , const gla_matrix * m1 , const gla_matrix * m2) ;
gla_matrix * gla_matrix_power(gla_matrix * m , const size_t k) ;
void gla_matrix_transpose(gla_matrix * m) ;

gla_vector * gla_matrix_rvector(const gla_matrix * m , gla_vector * v , const size_t i) ;
gla_vector * gla_matrix_cvector(const gla_matrix * m , gla_vector * v , const size_t j) ;

gla_matrix * gla_matrix_from_vector(gla_matrix * m , const gla_vector * pv) ;

int gla_matrix_is_invertible(const gla_matrix * m) ;
gla_matrix * gla_matrix_adjoint(const gla_matrix * m , gla_matrix * ad) ;
gla_matrix * gla_matrix_invertible(gla_matrix * des , const gla_matrix * src) ;

#endif // GLA_H
