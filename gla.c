#include "gla.h"
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

gla_determinant * gla_determinant_alloc(const size_t n) 
{
	gla_determinant * d ;
	if(0 == n) PrintErr("Determinant dimension should be a positive integer" , GLA_EINVAL) ;
	
	d = (gla_determinant *) malloc(sizeof(gla_determinant)) ;
	if(0 == d) PrintErr("Determinant alloc failure" , GLA_ENOMEM) ;
	
	d->dat = (type *) malloc(n * n * sizeof(type)) ;
	if(0 == d->dat) PrintErr("Determinant dat alloc failure" , GLA_ENOMEM) ;
	
	d->n = n ;
	return d ;
}

gla_determinant * gla_determinant_cpy(gla_determinant * c , const gla_determinant * d)
{
	if(0 == d || 0 == c) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	if(c->n != d->n) PrintErr("Dest gla_determinant 's size doesn't match src gla_determinant 's size" , GLA_EINVAL) ;
	
	memcpy(c->dat , d->dat , d->n * d->n * sizeof(type)) ;
	
	return c ;
}

void gla_determinant_set(gla_determinant * d , const size_t i , const size_t j , const type elem) 
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	if(d->n <= i || d->n <= j) PrintErr("Determinant doesn't have this elem" , GLA_EINVAL) ;
	
	d->dat[i * d->n + j] = elem ;
}

type gla_determinant_get(const gla_determinant * d , const size_t i , const size_t j) 
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	if(d->n <= i || d->n <= j) PrintErr("Determinant doesn't have this elem" , GLA_EINVAL) ;
	
	return d->dat[i * d->n + j] ;
}

void gla_determinant_free(gla_determinant * d) 
{
	if(d)
	{
		if(d->dat) free(d->dat) ;
		free(d) ;
	}
}

void gla_determinant_multiplication(gla_determinant * d , const type k , const size_t i)
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	if(i >= d->n) PrintErr("Row num should be smaller than gla_determinant dimension" , GLA_EINVAL) ;
	
	size_t j ;
	for(j = 0 ; j < d->n ; ++j)
		d->dat[i*d->n+j] *= k ;
}
void gla_determinant_addition(gla_determinant * d , const size_t des , const size_t src)
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	if(des >= d->n || src >= d->n) 
		PrintErr("Row num should be smaller than gla_determinant dimension" , GLA_EINVAL) ;
		
	size_t j ;
	type temp ;
	for(j = 0 ; j < d->n ; ++j)
	{
		temp = gla_determinant_get(d , des , j) + gla_determinant_get(d , src , j) ;
		gla_determinant_set(d , des , j , temp) ;
	}
}
void gla_determinant_swap(gla_determinant * d , const size_t i , const size_t i2)
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	if(i >= d->n || i2 >= d->n) 
		PrintErr("Row num should be smaller than gla_determinant dimension" , GLA_EINVAL) ;
	size_t j ;
	type temp ;
	for(j = 0 ; j < d->n ; ++j)
	{
		temp = gla_determinant_get(d , i , j) ;
		gla_determinant_set(d , i , j , gla_determinant_get(d , i2 , j)) ;
		gla_determinant_set(d , i2 , j , temp) ;
	}
}

type gla_determinant_alge_complement(const gla_determinant * d , const size_t i , const size_t j) 
{
	type result ;
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	if(d->n <= i || d->n <= j) PrintErr("Determinant doesn't have this elem" , GLA_EINVAL) ;
	
	if(1 == d->n) PrintErr("Determinant doesn't has alge complement" , GLA_EINVAL) ;
	
	gla_determinant * c = gla_determinant_alloc(d->n - 1) ;
	
	size_t x , y ;
	for(x = 0 , y = 0 ; x < d->n * d->n ; ++x)
	{
		if(i != x / d->n && j != x % d->n)
		{
            gla_determinant_set(c , y / c->n , y % c->n , gla_determinant_get(d , x / d->n , x % d->n)) ;
			++y ;
		}
	}
	assert(y == c->n * c->n) ;
	
	result = gla_determinant_value(c) ;
	
	if((i + j) % 2)
		result *= -1 ;
	
	gla_determinant_free(c) ;
	return result ;
}

type gla_determinant_value(const gla_determinant * d) 
{
	type result = 1 ;
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	gla_determinant * c = gla_determinant_alloc(d->n) ;
	
	c = gla_determinant_cpy(c , d) ;
	
	size_t i , j , k ;
	for(i = 0 ; i < c->n ; ++i)
	{
		type temp1 , temp2 , mul ;
		temp1 = gla_determinant_get(c , i , i) ;
		for(j = c->n - 1 ; j != i ; --j)
		{
			temp2 = gla_determinant_get(c , j , i) ;
			mul = -1 * temp2 / temp1 ;
			for(k = i ; k < c->n ; ++k)
			{
				gla_determinant_set(c , j , k , mul * gla_determinant_get(c , i , k) + gla_determinant_get(c , j , k)) ;
			}
		}
		
		for(j = i + 1 ; j < c->n ; ++j)
			assert(0 == gla_determinant_get(c , j , i)) ;
		result *= gla_determinant_get(c , i , i) ;
	}
	gla_determinant_free(c) ;
	return result ;
}

int gla_determinant_is_diagonal(const gla_determinant * d) 
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	return gla_determinant_is_upper_triangular(d) && gla_determinant_is_lower_triangular(d) ;
}
int gla_determinant_is_upper_triangular(const gla_determinant * d) 
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	size_t i , j ;
	for(i = 0 ; i < d->n ; ++i)
	{
		for(j = i + 1 ; j < d->n ; ++j)
		{
			if(0 != gla_determinant_get(d , i , j)) return 0 ;
		}
	}
	return 1 ;
}
int gla_determinant_is_lower_triangular(const gla_determinant * d) 
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	size_t i , j ;
	for(i = 0 ; i < d->n ; ++i)
	{
		for(j = 0 ; j < i ; ++j)
		{
			if(0 != gla_determinant_get(d , i , j)) return 0 ;
		}
	}
	return 1 ;
}
void gla_determinant_transpose(const gla_determinant * d) 
{
	if(0 == d) PrintErr("Determinant can't be null" , GLA_EINVAL) ;
	
	size_t i , j ;
	type temp ;
	for(i = 0 ; i < d->n ; ++i)
	{
		for(j = 0 ; j < i ; ++j)
		{
			temp = gla_determinant_get(d , i , j) ;
			gla_determinant_set(d , i , j , gla_determinant_get(d , j , i)) ;
			gla_determinant_set(d , j , i , temp) ;
		}
	}
}

//******************************************************

gla_vector * gla_vector_alloc(const size_t n) 
{
	gla_vector * result ;
	if(0 == n) PrintErr("Num should be a positive integer" , GLA_EINVAL) ;
	
	result = (gla_vector *) malloc(sizeof(gla_vector)) ;
	if(0 == result) PrintErr("Malloc failure" , GLA_EINVAL) ;
	
	result->dat = (type *) malloc(n * sizeof(type)) ;
	if(0 == result->dat) PrintErr("Malloc failure" , GLA_EINVAL) ;
	
	result->n = n ;
	
	return result ;
}

void gla_vector_set(gla_vector * v , const size_t j , type value) 
{
	if(0 == v) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	if(v->n <= j) PrintErr("Subscript should be smaller than Vector dimension" , GLA_EINVAL) ;
	
	v->dat[j] = value ;
}

type gla_vector_get(const gla_vector * v , const size_t j)
{
	if(0 == v) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	if(v->n <= j) PrintErr("Subscript should be smaller than Vector dimension" , GLA_EINVAL) ;
	
	return v->dat[j] ;
}
void gla_vector_free(gla_vector * v) 
{
	if(v)
	{
		if(v->dat) free(v->dat) ;
		free(v) ;
	}
}

gla_vector * gla_vector_cpy(gla_vector * des , const gla_vector * src) 
{
	if(0 == des || 0 == src) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	
	if(des->n != src->n) PrintErr("Des dimension's size isn't equal to src dimension's size" , GLA_EINVAL) ;
	
	memcpy(des->dat , src->dat , des->n * sizeof(type)) ;
	
	return des ;
}

void gla_vector_multiplication(gla_vector * v , const type k)
{
	if(0 == v) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	
	size_t j ;
	for(j = 0 ; j < v->n ; ++j)
		v->dat[j] *= k ;
}

void gla_vector_addition(gla_vector * des , const gla_vector * src)
{
	if(0 == des || 0 == src) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	if(des->n != src->n) PrintErr("Des dimension's size isn't equal to src dimension's size" , GLA_EINVAL) ;
	
	size_t j ;
	for(j = 0 ; j < des->n ; ++j)
		des->dat[j] += src->dat[j] ;
}
void gla_vector_swap(gla_vector * v1 , gla_vector * v2) 
{
	if(0 == v1 || 0 == v2) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	if(v1->n != v2->n) PrintErr("V1 dimension's size isn't equal to V2 dimension's size" , GLA_EINVAL) ;
	
	size_t j ;
	type temp ;
	for(j = 0 ; j < v1->n ; ++j)
	{
		temp = v1->dat[j] ;
		v1->dat[j] = v2->dat[j] ;
		v2->dat[j] = temp ;
	}
}
type gla_vector_dot_multiplication(const gla_vector * v1 , const gla_vector * v2)
{
	type result = 0 ;
	if(0 == v1 || 0 == v2) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	if(v1->n != v2->n) PrintErr("V1 dimension's size isn't equal to V2 dimension's size" , GLA_EINVAL) ;
	
	size_t j ;
	for(j = 0 ; j < v1->n ; ++j)
		result += v1->dat[j] * v2->dat[j] ;
	return result ;
}
gla_vector * gla_vector_cross_multiplication(const gla_vector * v1 , const gla_vector * v2 , gla_vector * des)
{
	if(0 == v1 || 0 == v2 || 0 == des) PrintErr("Vector can't be a null" , GLA_EINVAL) ;
	
	if(3 != v1->n || 3 != v2->n || 3 != des->n) PrintErr("Size isn't 3" , GLA_EINVAL) ;
	
	des->dat[0] = v1->dat[1] * v2->dat[2] - v1->dat[2] * v2->dat[1] ;
	des->dat[1] = -1 * (v1->dat[0] * v2->dat[2] - v1->dat[2] * v2->dat[0]) ;
	des->dat[2] = v1->dat[0] * v2->dat[1] - v1->dat[1] * v2->dat[0] ;
	return des ;
}

gla_matrix * gla_matrix_alloc(const size_t n , const size_t m)
{
    gla_matrix * result ;
    if(0 == n || 0 == m) PrintErr("Matrix dimension should be a positive integer" , GLA_EINVAL) ;
    
    result = (gla_matrix *) malloc(sizeof(gla_matrix)) ;
    if(0 == result) PrintErr("Malloc Failure" , GLA_ENOMEM) ;
    
    result->dat = (type *) malloc(n * m * sizeof(type)) ;
    if(0 == result->dat) PrintErr("Malloc Failure" , GLA_ENOMEM) ;
    
    result->n = n ;
    result->m = m ;
    return result ;
}
void gla_matrix_set(gla_matrix * m , const size_t i , const size_t j , const type ele)
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(i >= m->n || j >= m->m) PrintErr("Matrix doesn't have this elem" , GLA_EINVAL) ;
    
    m->dat[i*m->m+j] = ele ;
}
type gla_matrix_get(const gla_matrix * m , const size_t i , const size_t j)
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(i >= m->n || j >= m->m) PrintErr("Matrix doesn't have this elem" , GLA_EINVAL) ;

    return m->dat[i*m->m+j] ;
}
void gla_matrix_free(gla_matrix * m)
{
    if(m)
    {
        if(m->dat) free(m->dat) ;
        free(m) ;
    }
}
gla_matrix * gla_matrix_cpy(gla_matrix * des , const gla_matrix * src)
{
    if(0 == des || 0 == src) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(des->n != src->n || des->m != src->m) PrintErr("Matrix size doesn't match" , GLA_EINVAL) ;
    
    memcpy(des->dat , src->dat , src->n * src->m * sizeof(type)) ;
    return des ;
}
void gla_matrix_addition(gla_matrix * des , const gla_matrix * src) 
{
    if(0 == des || 0 == src) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(des->n != src->n || des->m != src->m) PrintErr("Matrix size doesn't match" , GLA_EINVAL) ;
    
    size_t i , j ;
    for(i = 0 ; i < des->n ; ++i)
    {
        for(j = 0 ; j < des->m ; ++j)
            des->dat[i*des->m+j] += src->dat[i*des->m+j] ;
    }
}
void gla_matrix_multiplication(gla_matrix * m , const type k) 
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    
    size_t i , j ;
    for(i = 0 ; i < m->n ; ++i)
    {
        for(j = 0 ; j < m->m ; ++j)
            m->dat[i*m->m+j] *= k ;
    }
}
gla_matrix * gla_matrix_matrix_multiplication(gla_matrix * des , const gla_matrix * m1 , const gla_matrix * m2) 
{
    if(0 == des || 0 == m1 || 0 == m2) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(m1->m != m2->n || m1->n != des->n || des->m != m2->m)
        PrintErr("Matrixs can't do multiplication operator" , GLA_EINVAL) ;
        
    size_t i , j , k ;
    for(i = 0 ; i < des->n ; ++i)
    {
        for(j = 0 ; j < des->m ; ++j)
        {
            des->dat[i*des->m+j] = 0 ;
            for(k = 0 ; k < m1->m ; ++k)
                des->dat[i*des->m+j] += m1->dat[i*m1->m+k] * m2->dat[k*m2->m+j] ;
        }
    }
    return des ;
}
gla_matrix * gla_matrix_square(gla_matrix * m)
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(m->n != m->m) PrintErr("Matrix 's row size != col size" , GLA_EINVAL) ;
    gla_matrix * mcopy = gla_matrix_alloc(m->n , m->m) ;
    
    gla_matrix_cpy(mcopy , m) ;
    
    m = gla_matrix_matrix_multiplication(m , mcopy , mcopy) ;
    
    gla_matrix_free(mcopy) ;
    return m ;
}
gla_matrix * gla_matrix_power(gla_matrix * m , const size_t k)
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(m->n != m->m) PrintErr("Matrix can't do this operator" , GLA_EINVAL) ;
    
    gla_matrix * temp = gla_matrix_alloc(m->n , m->m) ;
    temp = gla_matrix_cpy(temp , m) ;
    
    size_t i ;
    for(i = k ; i > 0 ; i >>= 1)
        gla_matrix_square(m) ;
    if(k % 2) 
    {
        gla_matrix * des = gla_matrix_alloc(m->n , m->m) ;
        des = gla_matrix_matrix_multiplication(des , m , temp) ;
        gla_matrix_cpy(m , des) ;
        gla_matrix_free(des) ;
    }
    gla_matrix_free(temp) ;
    return m ;
}
void gla_matrix_transpose(gla_matrix * m) 
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    
    size_t i , j ;
    type temp ;
    for(i = 0 ; i < m->n ; ++i)
    {
        for(j = 0 ; j < i ; ++j)
        {
            temp = m->dat[i*m->m+j] ;
            m->dat[i*m->m+j] = m->dat[j*m->m+i] ;
            m->dat[j*m->m+i] = temp ;
        }
    }
}

gla_vector * gla_matrix_rvector(const gla_matrix * m , gla_vector * v , const size_t i)
{
    if(0 == m || 0 == v) PrintErr("Matrix and Vector can't be a null" , GLA_EINVAL) ;
    if(m->m != v->n) PrintErr("Size doesn't match" , GLA_EINVAL) ;
    
    size_t j ;
    for(j = 0 ; j < m->m ; ++j)
        v->dat[j] = m->dat[i*m->m+j] ;
    return v ;
}
gla_vector * gla_matrix_cvector(const gla_matrix * m , gla_vector * v , const size_t j)
{
    if(0 == m || 0 == v) PrintErr("Matrix and Vector can't be a null" , GLA_EINVAL) ;
    if(m->m != v->n) PrintErr("Size doesn't match" , GLA_EINVAL) ;
    
    size_t i ;
    for(i = 0 ; i < m->n ; ++i)
        v->dat[i] = m->dat[i*m->m+j] ;
    return v ;
}

gla_matrix * gla_matrix_from_rvector(gla_matrix * m , const gla_vector * v)
{
    if(0 == m || 0 == v) PrintErr("Matrix and Vector can't be a null" , GLA_EINVAL) ;
    
    if(m->m != v->n) PrintErr("Size doesn't match" , GLA_EINVAL) ;
    
    size_t i , j ;
    for(i = 0 ; i < m->n ; ++i)
    {
        for(j = 0 ; j < m->m ; ++j)
        {
            m->dat[i*m->m+j] = v[i].dat[j] ;
        }
    }
    return m ;
}
gla_matrix * gla_matrix_from_cvector(gla_matrix * m , const gla_vector * v)
{
    if(0 == m || 0 == v) PrintErr("Matrix and Vector can't be a null" , GLA_EINVAL) ;
    if(m->n != v->n) PrintErr("Size doesn't match" , GLA_EINVAL) ;
    
    size_t i , j ;
    for(j = 0 ; j < m->m ; ++j)
    {
        for(i = 0 ; i < m->n ; ++i) 
            m->dat[i*m->m+j] = v[j].dat[i] ;
    }
    return m ;
}
int gla_matrix_is_invertible(const gla_matrix * m) 
{
    if(0 == m) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    
    if(m->n != m->m) PrintErr("Size doesn't match" , GLA_EINVAL) ;
    
    gla_determinant * d = gla_determinant_alloc(m->n) ;
    
    memcpy(d->dat , m->dat , m->n * m->n * sizeof(type)) ;
    
    type result = gla_determinant_value(d) ;
    
    gla_determinant_free(d) ;
    if(result > 0 || result < 0) return 1 ;
    else return 0 ;
}
gla_matrix * gla_matrix_adjoint(const gla_matrix * m , gla_matrix * ad) 
{
    if(0 == m || 0 == ad) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(m->n != m->m || ad->n != ad->m || m->n != ad->n) PrintErr("Size doesn't match" , GLA_EINVAL) ;
    
    gla_determinant * temp = gla_determinant_alloc(m->n) ;
    
    memcpy(temp->dat , m->dat , m->n * m->m * sizeof(type)) ;
    
    size_t i , j ;
    for(i = 0 ; i < m->n ; ++i)
    {
        for(j = 0 ; j < m->m ; ++j)
            ad->dat[i*ad->m+j] = gla_determinant_alge_complement(temp , i , j) ;
    }
    
    gla_determinant_free(temp) ;
    return ad ;
}
gla_matrix * gla_matrix_invertible(gla_matrix * des , const gla_matrix * src)
{
    if(0 == des || 0 == src) PrintErr("Matrix can't be a null" , GLA_EINVAL) ;
    if(des->n != des->m || src->n != src->m || des->n != src->n) 
        PrintErr("Size doesn't match" , GLA_EINVAL) ; 
    gla_matrix_adjoint(src , des) ;
    
    gla_determinant * d = gla_determinant_alloc(src->n) ;
    
    memcpy(d->dat , src->dat , src->n * src->n * sizeof(type)) ;
    
    type val = gla_determinant_value(d) ;
    
    gla_matrix_multiplication(des , 1 / val) ;
    
    gla_determinant_free(d) ;
    return des ;
}
