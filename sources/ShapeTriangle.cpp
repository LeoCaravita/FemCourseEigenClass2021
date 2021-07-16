//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"
#include "Shape1d.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1) {
        std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    auto nf = NShapeFunctions(orders);

    if (orders[nf-1] > 2) {
        std::cout << "ShapeTriangle::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }

    //número de funções
    int nshape = NShapeFunctions(orders);
    int nsides = nSides;

    //Define o tamanho de phi e dphi
    phi.resize(nshape);
    dphi.resize(2, nshape);

    // Linear order
    phi[0] = 1. - xi[0] - xi[1];
    phi[1] = xi[0];
    phi[2] = xi[1];

    dphi(0, 0) = -1.;
    dphi(1, 0) = -1.;
    dphi(0, 1) = 1.;
    dphi(1, 1) = 0.;
    dphi(0, 2) = 0.;
    dphi(1, 2) = 1.;

    int count = 3;

    //calcula os sides shapes
    for (int i = 3; i < 6; i++) {
        if (orders(i) == 2) {
            switch (count)
            {
            case 3:
                phi(count) = 4. * phi[0] * phi[1];
                dphi(0, count) = (1. - 2. * xi[0] * xi[1]) * 4.;
                dphi(1, count) = (xi[0] - xi[0] * xi[0]) * 4.;
                break;
            case 4:
                phi(count) = 4. * phi[1] * phi[2];
                dphi(0, count) = xi[1] * 4.;
                dphi(1, count) = xi[0] * 4.;
                break;
            case 5:
                phi(count) = 4. * phi[0] * phi[2];
                dphi(0, count) = (xi[1] - xi[1] * xi[1]) * 4.;
                dphi(1, count) = (-xi[1] * xi[1]) * 4.;
                break;
            }
            count++;
        }
        else if (orders(i) != 1) DebugStop();
    }


    //calcula o shape do internal side
    if (orders(6) == 3) {
        phi(count) = 27. * phi[0] * phi[1] * phi[2];
        dphi(0, count) = (xi[1] - 2. * xi[1] * xi[1] * xi[0]) * 27.;
        dphi(1, count) = (xi[0] - 2. * xi[0] * xi[0] * xi[1]) * 27.;
        count++;
    }
    else if (orders(6) != 1 && orders(6) != 2) DebugStop();


    if (count != nshape) DebugStop();
  
}

/// returns the number of shape functions associated with a side
int ShapeTriangle::NShapeFunctions(int side, int order){
    switch(side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
            return order-1;
        case 6:
            return 0;
    }
    
    DebugStop();
    std::cout << "ShapeTriangle::NShapeFunctions, bad parameter side " << std::endl;
    return 0;
}

/// returns the total number of shape functions
int ShapeTriangle::NShapeFunctions(VecInt &orders){
    
    int res=3;
    for(int in=3; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
    
}
