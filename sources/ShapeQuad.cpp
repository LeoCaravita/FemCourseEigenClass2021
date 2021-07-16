//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    auto nf = NShapeFunctions(orders);

    if (orders[nf-1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }    

    //número de funções
    int nshape = NShapeFunctions(orders);
    int nsides = nSides;

    //Define o tamanho de phi e dphi
    phi.resize(nshape);
    dphi.resize(2, nshape);

    //Calcula os phis lineares
    phi(0) = (1. - xi[0]) * (1. - xi[1]) * (1. / 4.);
    phi(1) = (1. + xi[0]) * (1. - xi[1]) * (1. / 4.);
    phi(2) = (1. + xi[0]) * (1. + xi[1]) * (1. / 4.);
    phi(3) = (1. - xi[0]) * (1. + xi[1]) * (1. / 4.);

    // calculas os dphis lineares
    dphi(0, 0) = (-1. + xi[1]) * (1. / 4.);
    dphi(0, 1) = (1. - xi[1]) * (1. / 4.);
    dphi(0, 2) = (1. + xi[1]) * (1. / 4.);
    dphi(0, 3) = (-1. - xi[1]) * (1. / 4.);

    dphi(1, 0) = (-1. + xi[0]) * (1. / 4.);
    dphi(1, 1) = (-1. - xi[0]) * (1. / 4.);
    dphi(1, 2) = (1. + xi[0]) * (1. / 4.);
    dphi(1, 3) = (1. - xi[0]) * (1. / 4.);

    int count = 4;

    //calcula os sides shapes
    for (int i = 4; i < 8; i++) {
        if (orders(i) == 2) {
            switch (count)
            {
            case 4:
                phi(count) = (1. - xi[0] * xi[0]) * (1. - xi[1]) / 2.;
                dphi(0, count) = (-2. * xi[0] + 2. * xi[0] * xi[1]) / 2.;
                dphi(1, count) = (-1. + xi[0] * xi[0]) / 2.;
                break;
            case 5:
                phi(count) = (1. + xi[0]) * (1. - xi[1] * xi[1]) / 2.;
                dphi(0, count) = (1. - xi[1] * xi[1]) / 2.;
                dphi(1, count) = (-2. * xi[1] - 2. * xi[1] * xi[0]) / 2.;
                break;
            case 6:
                phi(count) = (1. - xi[0] * xi[0]) * (1. + xi[1]) / 2.;
                dphi(0, count) = (-2. * xi[0] - 2. * xi[0] * xi[1]) / 2.;
                dphi(1, count) = (1. - xi[0] * xi[0]) / 2.;
                break;
            case 7:
                phi(count) = (1. - xi[0]) * (1. - xi[1] * xi[1]) / 2.;
                dphi(0, count) = (-1. + xi[1] * xi[1]) / 2.;
                dphi(1, count) = (-2. * xi[1] + 2. * xi[1] * xi[0]) / 2.;
                break;
            }
            count++;
        }
        else if (orders(i) != 1) DebugStop();
    }


    //calcula o shape do internal side
    if (orders(8) == 2) {
        phi(count) = (1. - xi[0] * xi[0]) * (1. - xi[1] * xi[1]);
        dphi(0, count) = (-2. * xi[0] + 2. * xi[1] * xi[1] * xi[0]);
        dphi(1, count) = (-2. * xi[1] + 2. * xi[1] * xi[0] * xi[0]);
        count++;
    }
    else if (orders(8) != 1) DebugStop();


    if (count != nshape) DebugStop();

}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
   
    return res;
}
