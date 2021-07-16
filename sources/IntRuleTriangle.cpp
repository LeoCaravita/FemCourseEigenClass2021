/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRuleTriangle.h"

IntRuleTriangle::IntRuleTriangle(){

}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    fOrder = order;

#define Perm3(a)	a, a, a
#define Dup3(w)		(0.5*w)
#define Perm21(a)	a, a, 1. - 2.0*a, a, 1.- 2.*a, a, 1.- 2.*a, a, a
#define Dup21(w)	Dup3(w), Dup3(w), Dup3(w)
#define Perm111(a,b)	a, b, 1.L - a - b, a, 1.L - a - b, b, b, a , 1.L - a - b, \
	    b, 1.L - a - b, a, 1.L - a - b, a, b, 1.L - a - b, b, a
#define Dup111(w)	Dup3(w), Dup3(w), Dup3(w), Dup3(w), Dup3(w), Dup3(w)

double pontos1[] = { Perm3(1. / 3) };

double pontos2[] = { Perm21(1. / 6) };
double pesos2[] = { Dup21(1. / 3) };
int contador2 = 0;

double pontos4[] = {
    Perm21(.44594849091596488631832925388305199),
    Perm21(.09157621350977074345957146340220151) };
double pesos4[] = {
    Dup21(.22338158967801146569500700843312280),
    Dup21(.10995174365532186763832632490021053) };
int contador4 = 0;

double pontos5[] = {
    Perm21(.10128650732345633880098736191512383),
    Perm21(.47014206410511508977044120951344760),
    Perm3(1. / 3) };
double pesos5[] = {
    Dup21(.12593918054482715259568394550018133),
    Dup21(.13239415278850618073764938783315200),
    Dup3(9. / 40) };
int contador5 = 0;


if (order < 0 || order > MaxOrder()) DebugStop();
int i;
int j;

switch (order)
{
case 0:
case 1:
    fPoints.resize(Dimension(), 3);
    fWeights.resize(1);

    for (i = 0; i <= 2;i++) {
        fPoints(0, i) = pontos1[i];
    }
    fWeights[0] = Dup3(1);
    break;

case 2:
    fPoints.resize(3, 3);
    fWeights.resize(3);
    for (i = 0; i <= 2;i++) {
        for (j = 0; j <= 2; j++) {
            fPoints(i, j) = pontos2[contador2];
            contador2++;
        }
    }
    for (i = 0; i <= 2;i++) {
        fWeights[i] = pesos2[i];
    }
    break;

case 3:
case 4:
    fPoints.resize(6, 3);
    fWeights.resize(6);
    for (i = 0; i <= 5;i++) {
        for (j = 0; j <= 2; j++) {
            fPoints(i, j) = pontos4[contador4];
            contador4++;
        }
    }
    for (i = 0; i <= 5;i++) {
        fWeights[i] = pesos4[i];
    }
    break;

case 5:
    fPoints.resize(7, 3);
    fWeights.resize(7);
    for (i = 0; i <= 6;i++) {
        for (j = 0; j <= 2; j++) {
            fPoints(i, j) = pontos5[contador5];
            contador5++;
        }
    }
    for (i = 0; i <= 6;i++) {
        fWeights[i] = pesos5[i];
    }
    break;
default:
    DebugStop();
    break;
}
}