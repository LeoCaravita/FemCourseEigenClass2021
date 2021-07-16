
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}


void Geom1d::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    if (xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) {
        DebugStop();
    }

    phi.resize(2);
    phi(0) = (1. - xi[0]) / 2.;
    phi(1) = (1. + xi[0]) / 2.;

    dphi.resize(Dimension, 2);
    dphi(0, 0) = -1. / 2.;
    dphi(0, 1) = 1. / 2.;
}

void Geom1d::X(const VecDouble& xi, MatrixDouble& NodeCo, VecDouble& x) {
    if (xi.size() != Dimension) {
        DebugStop();
    }
    if (x.size() < NodeCo.rows()) {
        DebugStop();
    }
    if (NodeCo.cols() != nCorners) {
        DebugStop();
    }

    //Chama os valores da Shape
    VecDouble phi(2);
    MatrixDouble dphi(Dimension, 2);
    Shape(xi, phi, dphi);

    //Coleta o número de linhas da matriz de coordenadas
    int number_of_rows = NodeCo.rows();

    //Calcula o valor de X para cada coordenada
    for (int i = 0; i < number_of_rows; i++) {
        x[i] = NodeCo(i, 0) * phi[0] + NodeCo(i, 1) * phi[1];
    }

}

void Geom1d::GradX(const VecDouble& xi, MatrixDouble& NodeCo, VecDouble& x, MatrixDouble& gradx) {
    if (xi.size() != Dimension) {
        DebugStop();
    }
    if (x.size() != NodeCo.rows()) {
        DebugStop();
    }
    if (NodeCo.cols() != nCorners) {
        DebugStop();
    }

    //Chama os valores da Shape
    VecDouble phi(2);
    MatrixDouble dphi(Dimension, 2);
    Shape(xi, phi, dphi);

    //Coleta o número de linhas e colunas da matriz de coordenadas
    int number_of_rows = NodeCo.rows();
    int number_of_columns = NodeCo.cols();

    //Define a dimensão do gradiente e de X, setando com zeros
    gradx.resize(number_of_rows, 1);
    gradx.setZero();
    x.resize(number_of_rows);
    x.setZero();

    //Calcula o valor do gradiente para cada coordenada
    for (int i = 0; i < number_of_columns; i++) {
        for (int j = 0; j < number_of_rows; j++) {
            x(j) += NodeCo(j, i) * phi(i);
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
        }
    }
}


void Geom1d::SetNodes(const VecInt& nodes) {

    if (nodes.rows() != 2) DebugStop();

    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt& nodes) const {
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) const {
    if (node < 0 || node > 2) DebugStop();
    return fNodeIndices[node];
}

int Geom1d::NumNodes() {
    return nCorners;
}

GeoElementSide Geom1d::Neighbour(int side) const {
    if (side < 0 || side>2) DebugStop();
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide& neighbour) {
    if (side < 0 || side > 2) DebugStop();
    fNeighbours[side] = neighbour;
}

