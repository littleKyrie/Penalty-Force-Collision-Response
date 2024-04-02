#include <d2psi_neo_hookean_dq2.h>
#include <iostream>

void d2psi_neo_hookean_dF2(
        Eigen::Matrix99d &dP,
        Eigen::Ref<const Eigen::Matrix3d> F,
        double C, double D) {

    Eigen::Matrix99d ddw;

    // prepare data
    double dxdX = F(0, 0);
    double dxdY = F(0, 1);
    double dxdZ = F(0, 2);

    double dydX = F(1, 0);
    double dydY = F(1, 1);
    double dydZ = F(1, 2);

    double dzdX = F(2, 0);
    double dzdY = F(2, 1);
    double dzdZ = F(2, 2);

    // assemble dP
    // the row index represent the index of F's argument(row first in F)
    // the col index represent the index of F's dependent variable(row first in P)
    // dP/dF00
    dP(0, 0) = C * (pow((dydY * dzdZ - dydZ * dzdY), 2) /
                    pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                         dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
             + (D * pow((dydY * dzdZ - dydZ * dzdY), 2)) /
               pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                    dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
             - (D * log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX) *
                    pow((dydY * dzdZ - dydZ * dzdY), 2)) /
               pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                    dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 0) = (D * log(dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX)
            * (dydX * dzdZ - dydZ * dzdX) * (dydY * dzdZ - dydZ * dzdY)) /
            pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
            - (D*(dydX*dzdZ - dydZ*dzdX)*(dydY*dzdZ - dydZ*dzdY)) /
            pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
            - (C*(dydX*dzdZ - dydZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))/
            pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX),2);
    dP(2,0) = (C * (dydX * dzdY - dydY * dzdX) * (dydY * dzdZ - dydZ * dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
            + (D * (dydX * dzdY - dydY * dzdX) * (dydY * dzdZ - dydZ * dzdY)) /
              pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
            - (D * log(dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX)
            * (dydX * dzdY - dydY * dzdX) * (dydY * dzdZ - dydZ * dzdY))/
              pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3,0) = (D * log(dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdY*dzdZ - dxdZ*dzdY)*(dydY*dzdZ - dydZ*dzdY))/
              pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - (D*(dxdY*dzdZ - dxdZ*dzdY)*(dydY*dzdZ - dydZ*dzdY))/
                     pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (C*(dxdY*dzdZ - dxdZ*dzdY)*(dydY*dzdZ - dydZ*dzdY))/
                            pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 0) = (D*(dxdX*dzdZ - dxdZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dzdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dzdZ - dxdZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dzdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dzdZ - dxdZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 0) = C*(dzdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dzdY - dxdY*dzdX)*(dydY*dzdZ - dydZ*dzdY))/
              pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dzdY - dxdY*dzdX)*(dydY*dzdZ - dydZ*dzdY))/
                             pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dzdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dzdY - dxdY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 0) = (C*(dxdY*dydZ - dxdZ*dydY)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdY*dydZ - dxdZ*dydY)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dydZ - dxdZ*dydY)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 0) = C*(dydZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dydZ - dxdZ*dydX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dydZ - dxdZ*dydX)*(dydY*dzdZ - dydZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dydZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydZ - dxdZ*dydX)*(dydY*dzdZ - dydZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 0) = (D*(dxdX*dydY - dxdY*dydX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dydY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dydY - dxdY*dydX)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dydY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    // dP/dF01
    dP(0, 1) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dydX*dzdZ - dydZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dydX*dzdZ - dydZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dydX*dzdZ - dydZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 1) = C*(pow((dydX * dzdZ - dydZ * dzdX), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dydX * dzdZ - dydZ * dzdX), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)*
                                       pow((dydX * dzdZ - dydZ * dzdX), 2))
                                       / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                              dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 1) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dydX*dzdY - dydY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dydX*dzdY - dydY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dydX*dzdY - dydY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 1) = C*(dzdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dzdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 1) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 1) = (D*(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dzdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dzdY - dxdY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dzdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 1) = (D*dydZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dydZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdY*dydZ - dxdZ*dydY)*(dydX*dzdZ - dydZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdZ - dydZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 1) = (C*(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdZ - dydZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdZ - dydZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 1) = C*(dydX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dydY - dxdY*dydX)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dydY - dxdY*dydX)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dydX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydY - dxdY*dydX)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    // dPdF02
    dP(0, 2) = (C*(dydX*dzdY - dydY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dydX*dzdY - dydY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dydX*dzdY - dydY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 2) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dydX*dzdY - dydY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dydX*dzdY - dydY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dydX*dzdY - dydY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 2) = C*(pow((dydX * dzdY - dydY * dzdX), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dydX * dzdY - dydY * dzdX), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)*
                                       pow((dydX * dzdY - dydY * dzdX), 2))
                                       / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                              dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 2) = (D*dzdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dzdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdY - dydY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdY - dydY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 2) = C*(dzdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdY - dydY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dzdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdY - dydY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 2) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdY - dydY*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdY - dydY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdY - dydY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 2) = C*(dydY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdY*dydZ - dxdZ*dydY)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdY - dydY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dydY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdY - dydY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 2) = (D*dydX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dydX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdX*dydZ - dxdZ*dydX)*(dydX*dzdY - dydY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdY - dydY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 2) = (C*(dxdX*dydY - dxdY*dydX)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dydY - dxdY*dydX)*(dydX*dzdY - dydY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dydX*dzdY - dydY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    //dP/dF10
    dP(0, 3) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdY*dzdZ - dxdZ*dzdY)*(dydY*dzdZ - dydZ*dzdY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdY*dzdZ - dxdZ*dzdY)*(dydY*dzdZ - dydZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdY*dzdZ - dxdZ*dzdY)*(dydY*dzdZ - dydZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 3) = C*(dzdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dzdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 3) = (D*dzdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dzdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdY - dydY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dzdZ - dxdZ*dzdY)*(dydX*dzdY - dydY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 3) = C*(pow((dxdY * dzdZ - dxdZ * dzdY), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dxdY * dzdZ - dxdZ * dzdY), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  * pow((dxdY * dzdZ - dxdZ * dzdY), 2))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 3) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdZ - dxdZ*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdZ - dxdZ*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdZ - dxdZ*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 3) = (C*(dxdX*dzdY - dxdY*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dzdY - dxdY*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dzdY - dxdY*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 3) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdY*dydZ - dxdZ*dydY)*(dxdY*dzdZ - dxdZ*dzdY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdY*dydZ - dxdZ*dydY)*(dxdY*dzdZ - dxdZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdY*dydZ - dxdZ*dydY)*(dxdY*dzdZ - dxdZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 3) = (D*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dxdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dydZ - dxdZ*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dxdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydZ - dxdZ*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 3) = C*(dxdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dydY - dxdY*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dydY - dxdY*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dxdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydY - dxdY*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    // dP/dF11
    dP(0, 4) = (D*(dxdX*dzdZ - dxdZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dzdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dzdZ - dxdZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dzdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dzdZ - dxdZ*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 4) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 4) = C*(dzdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdY - dydY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) -
                                  (D*dzdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dzdZ - dxdZ*dzdX)*(dydX*dzdY - dydY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 4) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdZ - dxdZ*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdZ - dxdZ*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdZ - dxdZ*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 4) = C*(pow((dxdX * dzdZ - dxdZ * dzdX), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dxdX * dzdZ - dxdZ * dzdX), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  * pow((dxdX * dzdZ - dxdZ * dzdX), 2))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 4) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdY - dxdY*dzdX)*(dxdX*dzdZ - dxdZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdY - dxdY*dzdX)*(dxdX*dzdZ - dxdZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdY - dxdY*dzdX)*(dxdX*dzdZ - dxdZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 4) = C*(dxdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdZ - dxdZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdZ - dxdZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dxdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdZ - dxdZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 4) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 4) = (D*(dxdX*dydY - dxdY*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dxdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dydY - dxdY*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dxdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    //dP/dF12
    dP(0, 5) = C*(dzdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dzdY - dxdY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dzdY - dxdY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dzdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dzdY - dxdY*dzdX)*(dydY*dzdZ - dydZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 5) = (D*(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dzdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dzdY - dxdY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dzdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdZ - dydZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 5) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdY - dydY*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdY - dydY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdY - dxdY*dzdX)*(dydX*dzdY - dydY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 5) = (C*(dxdX*dzdY - dxdY*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dzdY - dxdY*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dzdY - dxdY*dzdX)*(dxdY*dzdZ - dxdZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 5) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dzdY - dxdY*dzdX)*(dxdX*dzdZ - dxdZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dzdY - dxdY*dzdX)*(dxdX*dzdZ - dxdZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dzdY - dxdY*dzdX)*(dxdX*dzdZ - dxdZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 5) = C*(pow((dxdX * dzdY - dxdY * dzdX), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dxdX * dzdY - dxdY * dzdX), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  * pow((dxdX * dzdY - dxdY * dzdX), 2))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 5) = (D*dxdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdY - dxdY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dxdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdY - dxdY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdY - dxdY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 5) = C*(dxdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdY - dxdY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdY - dxdY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dxdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdY - dxdY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 5) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydY - dxdY*dydX)*(dxdX*dzdY - dxdY*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydY - dxdY*dydX)*(dxdX*dzdY - dxdY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydY - dxdY*dydX)*(dxdX*dzdY - dxdY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    // dP/dF20
    dP(0, 6) = (C*(dxdY*dydZ - dxdZ*dydY)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdY*dydZ - dxdZ*dydY)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dydZ - dxdZ*dydY)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 6) = (D*dydZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dydZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdY*dydZ - dxdZ*dydY)*(dydX*dzdZ - dydZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdZ - dydZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 6) = C*(dydY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdY*dydZ - dxdZ*dydY)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdY - dydY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dydY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdY*dydZ - dxdZ*dydY)*(dydX*dzdY - dydY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 6) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdY*dydZ - dxdZ*dydY)*(dxdY*dzdZ - dxdZ*dzdY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdY*dydZ - dxdZ*dydY)*(dxdY*dzdZ - dxdZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdY*dydZ - dxdZ*dydY)*(dxdY*dzdZ - dxdZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 6) = C*(dxdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdZ - dxdZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdZ - dxdZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dxdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdZ - dxdZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 6) = (D*dxdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdY - dxdY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dxdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdY - dxdY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdY*dydZ - dxdZ*dydY)*(dxdX*dzdY - dxdY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 6) = C*(pow((dxdY * dydZ - dxdZ * dydY), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dxdY * dydZ - dxdZ * dydY), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  * pow((dxdY * dydZ - dxdZ * dydY), 2))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 6) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dydZ - dxdZ*dydY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - (D*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dydZ - dxdZ*dydY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (C*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dydZ - dxdZ*dydY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 6) = (C*(dxdX*dydY - dxdY*dydX)*(dxdY*dydZ - dxdZ*dydY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dydY - dxdY*dydX)*(dxdY*dydZ - dxdZ*dydY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dxdY*dydZ - dxdZ*dydY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    // dP/dF21
    dP(0, 7) = C*(dydZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dydZ - dxdZ*dydX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dydZ - dxdZ*dydX)*(dydY*dzdZ - dydZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dydZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydZ - dxdZ*dydX)*(dydY*dzdZ - dydZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 7) = (C*(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdZ - dydZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdZ - dydZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 7) = (D*dydX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
            /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - (D*(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dydX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   + ((dxdX*dydZ - dxdZ*dydX)*(dydX*dzdY - dydY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydZ - dxdZ*dydX)*(dydX*dzdY - dydY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 7) = (D*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dxdZ/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dydZ - dxdZ*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dxdZ*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydZ - dxdZ*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 7) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 7) = C*(dxdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            + ((dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdY - dxdY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           + (D*(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdY - dxdY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dxdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX) -
                                  (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydZ - dxdZ*dydX)*(dxdX*dzdY - dxdY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 7) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydZ - dxdZ*dydX)*(dxdY*dydZ - dxdZ*dydY))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dydZ - dxdZ*dydY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydZ - dxdZ*dydX)*(dxdY*dydZ - dxdZ*dydY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 7) = C*(pow((dxdX * dydZ - dxdZ * dydX), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dxdX * dydZ - dxdZ * dydX), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  * pow((dxdX * dydZ - dxdZ * dydX), 2))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 7) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydY - dxdY*dydX)*(dxdX*dydZ - dxdZ*dydX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydY - dxdY*dydX)*(dxdX*dydZ - dxdZ*dydX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydY - dxdY*dydX)*(dxdX*dydZ - dxdZ*dydX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);

    // dPdF22
    dP(0, 8) = (D*(dxdX*dydY - dxdY*dydX)*(dydY*dzdZ - dydZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dydY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dydY - dxdY*dydX)*(dydY*dzdZ - dydZ*dzdY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dydY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dydY*dzdZ - dydZ*dzdY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(1, 8) = C*(dydX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dydY - dxdY*dydX)*(dydX*dzdZ - dydZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dydY - dxdY*dydX)*(dydX*dzdZ - dydZ*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dydX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydY - dxdY*dydX)*(dydX*dzdZ - dydZ*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(2, 8) = (C*(dxdX*dydY - dxdY*dydX)*(dydX*dzdY - dydY*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dydY - dxdY*dydX)*(dydX*dzdY - dydY*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dydX*dzdY - dydY*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(3, 8) = C*(dxdY/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            - ((dxdX*dydY - dxdY*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                           - (D*(dxdX*dydY - dxdY*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*dxdY*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                                  /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  + (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  *(dxdX*dydY - dxdY*dydX)*(dxdY*dzdZ - dxdZ*dzdY))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(4, 8) = (D*(dxdX*dydY - dxdY*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   - C*(dxdX/(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                   - ((dxdX*dydY - dxdY*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2))
                          + (D*dxdX*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX))
                          /(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dxdX*dzdZ - dxdZ*dzdX))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(5, 8) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydY - dxdY*dydX)*(dxdX*dzdY - dxdY*dzdX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydY - dxdY*dydX)*(dxdX*dzdY - dxdY*dzdX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydY - dxdY*dydX)*(dxdX*dzdY - dxdY*dzdX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(6, 8) = (C*(dxdX*dydY - dxdY*dydX)*(dxdY*dydZ - dxdZ*dydY))
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                   + (D*(dxdX*dydY - dxdY*dydX)*(dxdY*dydZ - dxdZ*dydY))
                   / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                          dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                          - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                          *(dxdX*dydY - dxdY*dydX)*(dxdY*dydZ - dxdZ*dydY))
                          / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                 dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(7, 8) = (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
            *(dxdX*dydY - dxdY*dydX)*(dxdX*dydZ - dxdZ*dydX))
                    / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                           dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                           - (D*(dxdX*dydY - dxdY*dydX)*(dxdX*dydZ - dxdZ*dydX))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (C*(dxdX*dydY - dxdY*dydX)*(dxdX*dydZ - dxdZ*dydX))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);
    dP(8, 8) = C*(pow((dxdX * dydY - dxdY * dydX), 2)
            / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                   dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2) + 1)
                           + (D* pow((dxdX * dydY - dxdY * dydX), 2))
                           / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ + dxdY * dydZ * dzdX +
                                  dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2)
                                  - (D*log(dxdX*dydY*dzdZ - dxdX*dydZ*dzdY - dxdY*dydX*dzdZ + dxdY*dydZ*dzdX + dxdZ*dydX*dzdY - dxdZ*dydY*dzdX)
                                  * pow((dxdX * dydY - dxdY * dydX), 2))
                                  / pow((dxdX * dydY * dzdZ - dxdX * dydZ * dzdY - dxdY * dydX * dzdZ +
                                         dxdY * dydZ * dzdX + dxdZ * dydX * dzdY - dxdZ * dydY * dzdX), 2);


    /*
    //-------------
    ddw(0,0) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F2_2 * F3_3 - F2_3 * F3_2, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) - F1_1 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F2_2 * F3_3 - F2_3 * F3_2, 2.0) * 2.0;
    ddw(0,1) = -C * (F1_1 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F1_2 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F2_1 * F3_3 - F2_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F2_1 * F3_3 - F2_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(0,2) = -C * (F1_1 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F1_3 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F2_1 * F3_2 - F2_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F2_1 * F3_2 - F2_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(0,3) = -C * (F1_1 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F2_1 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F3_3 - F1_3 * F3_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F3_3 - F1_3 * F3_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(0,4) = -C * (F3_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_2 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F3_3 - F1_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_3 - F1_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 - D * F3_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(0,5) = C * (F3_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_3 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F3_2 - F1_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_2 - F1_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 + D * F3_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(0,6) = -C * (F1_1 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_1 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_2 * F2_3 - F1_3 * F2_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F2_3 - F1_3 * F2_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(0,7) = C * (F2_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_2 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_3 - F1_3 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 + D * F2_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(0,8) = -C * (F2_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 - D * F2_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(1,0) = -C * (F1_1 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F1_2 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F2_1 * F3_3 - F2_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F2_1 * F3_3 - F2_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(1,1) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F2_1 * F3_3 - F2_3 * F3_1, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) + F1_2 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F2_1 * F3_3 - F2_3 * F3_1, 2.0) * 2.0;
    ddw(1,2) = -C * (F1_2 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F1_3 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F2_1 * F3_2 - F2_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F2_1 * F3_2 - F2_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0;
    ddw(1,3) = C * (F3_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_2 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_1 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 + D * F3_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(1,4) = -C * (F1_2 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_2 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0;
    ddw(1,5) = C * (F3_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (-2.0 / 3.0) + F1_2 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_3 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 - D * F3_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(1,6) = -C * (F2_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_2 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 - D * F2_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(1,7) = C * (F1_2 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0;
    ddw(1,8) = C * (F2_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_2 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 + D * F2_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(2,0) = -C * (F1_1 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F1_3 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F2_1 * F3_2 - F2_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F2_1 * F3_2 - F2_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(2,1) = -C * (F1_2 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F1_3 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F2_1 * F3_2 - F2_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F2_1 * F3_2 - F2_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0;
    ddw(2,2) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F2_1 * F3_2 - F2_2 * F3_1, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) - F1_3 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F2_1 * F3_2 - F2_2 * F3_1, 2.0) * 2.0;
    ddw(2,3) = -C * (F3_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_1 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 - D * F3_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(2,4) = C * (F3_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_2 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 + D * F3_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(2,5) = -C * (F1_3 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F2_3 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0;
    ddw(2,6) = C * (F2_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 + D * F2_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(2,7) = -C * (F2_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 - D * F2_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(2,8) = -C * (F1_3 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0;
    ddw(3,0) = -C * (F1_1 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F2_1 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F3_3 - F1_3 * F3_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F3_3 - F1_3 * F3_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(3,1) = C * (F3_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_2 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_1 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 + D * F3_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(3,2) = -C * (F3_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_1 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F3_3 - F1_3 * F3_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 - D * F3_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(3,3) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F1_2 * F3_3 - F1_3 * F3_2, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) + F2_1 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F1_2 * F3_3 - F1_3 * F3_2, 2.0) * 2.0;
    ddw(3,4) = -C * (F2_1 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_2 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_3 - F1_3 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_3 - F1_3 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0;
    ddw(3,5) = C * (F2_1 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_3 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_2 - F1_2 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0;
    ddw(3,6) = -C * (F2_1 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0;
    ddw(3,7) = C * (F1_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (-2.0 / 3.0) + F2_1 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0 - D * F1_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(3,8) = C * (F1_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F2_1 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0 + D * F1_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(4,0) = -C * (F3_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_2 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F3_3 - F1_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_3 - F1_3 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 - D * F3_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(4,1) = -C * (F1_2 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_2 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0;
    ddw(4,2) = C * (F3_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_2 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_3 - F1_3 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 + D * F3_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(4,3) = -C * (F2_1 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_2 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_3 - F1_3 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_3 - F1_3 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0;
    ddw(4,4) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F1_1 * F3_3 - F1_3 * F3_1, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) - F2_2 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F1_1 * F3_3 - F1_3 * F3_1, 2.0) * 2.0;
    ddw(4,5) = -C * (F2_2 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F2_3 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_2 - F1_2 * F3_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0;
    ddw(4,6) = C * (F1_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F2_2 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0 + D * F1_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(4,7) = -C * (F2_2 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F3_2 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0;
    ddw(4,8) = -C * (F1_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F2_2 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0 - D * F1_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(5,0) = C * (F3_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F2_3 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F3_2 - F1_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_2 - F1_2 * F3_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 + D * F3_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(5,1) = C * (F3_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (-2.0 / 3.0) + F1_2 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_3 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 - D * F3_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(5,2) = -C * (F1_3 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F2_3 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_2 - F1_2 * F3_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0;
    ddw(5,3) = C * (F2_1 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F2_3 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F3_2 - F1_2 * F3_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0;
    ddw(5,4) = -C * (F2_2 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F2_3 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F3_2 - F1_2 * F3_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F3_2 - F1_2 * F3_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0;
    ddw(5,5) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F1_1 * F3_2 - F1_2 * F3_1, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) + F2_3 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F1_1 * F3_2 - F1_2 * F3_1, 2.0) * 2.0;
    ddw(5,6) = -C * (F1_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F2_3 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_2 - F1_2 * F3_1) * 2.0 - D * F1_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(5,7) = C * (F1_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F2_3 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 2.0 + D * F1_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(5,8) = -C * (F2_3 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_3 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 2.0;
    ddw(6,0) = -C * (F1_1 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_1 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_2 * F2_3 - F1_3 * F2_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F2_3 - F1_3 * F2_2) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0;
    ddw(6,1) = -C * (F2_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_2 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 - D * F2_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(6,2) = C * (F2_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F2_3 - F1_3 * F2_2) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 + D * F2_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(6,3) = -C * (F2_1 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0;
    ddw(6,4) = C * (F1_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F2_2 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0 + D * F1_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(6,5) = -C * (F1_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F2_3 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_1 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * (F1_1 * F3_2 - F1_2 * F3_1) * 2.0 - D * F1_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(6,6) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F1_2 * F2_3 - F1_3 * F2_2, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) - F3_1 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F1_2 * F2_3 - F1_3 * F2_2, 2.0) * 2.0;
    ddw(6,7) = -C * (F3_1 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F3_2 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 2.0;
    ddw(6,8) = -C * (F3_1 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 2.0;
    ddw(7,0) = C * (F2_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_2 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_3 - F1_3 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 + D * F2_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(7,1) = C * (F1_2 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0;
    ddw(7,2) = -C * (F2_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_3 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0 - D * F2_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(7,3) = C * (F1_3 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (-2.0 / 3.0) + F2_1 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0 - D * F1_3 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(7,4) = -C * (F2_2 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F3_2 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0;
    ddw(7,5) = C * (F1_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F2_3 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_2 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 2.0 + D * F1_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(7,6) = -C * (F3_1 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (-4.0 / 3.0) + F3_2 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_3 - F1_3 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 2.0;
    ddw(7,7) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F1_1 * F2_3 - F1_3 * F2_1, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) + F3_2 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F1_1 * F2_3 - F1_3 * F2_1, 2.0) * 2.0;
    ddw(7,8) = -C * (F3_2 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_3 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F2_3 - F1_3 * F2_1) * 2.0;
    ddw(8,0) = -C * (F2_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F1_1 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F2_2 * F3_3 - F2_3 * F3_2) * 2.0 - D * F2_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(8,1) = C * (F2_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F1_2 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_3 - F2_3 * F3_1) * 2.0 + D * F2_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(8,2) = -C * (F1_3 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F2_1 * F3_2 - F2_2 * F3_1) * 2.0;
    ddw(8,3) = C * (F1_2 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) - F2_1 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F3_3 - F1_3 * F3_2) * 2.0 + D * F1_2 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(8,4) = -C * (F1_1 * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (2.0 / 3.0) + F2_2 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_3 - F1_3 * F3_1) * 2.0 - D * F1_1 * (-F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0) * 2.0;
    ddw(8,5) = -C * (F2_3 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_3 * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F3_2 - F1_2 * F3_1) * 2.0;
    ddw(8,6) = -C * (F3_1 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + F3_3 * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) + D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_2 * F2_3 - F1_3 * F2_2) * 2.0;
    ddw(8,7) = -C * (F3_2 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) - F3_3 * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (4.0 / 3.0) + (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * (F1_1 * F2_3 - F1_3 * F2_1) * 2.0;
    ddw(8,8) = C * (1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 2.0 / 3.0) * 2.0 + pow(F1_1 * F2_2 - F1_2 * F2_1, 2.0) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 8.0 / 3.0) * (F1_1 * F1_1 + F1_2 * F1_2 + F1_3 * F1_3 + F2_1 * F2_1 + F2_2 * F2_2 + F2_3 * F2_3 + F3_1 * F3_1 + F3_2 * F3_2 + F3_3 * F3_3) * (1.0E1 / 9.0) - F3_3 * (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F1_1 * F2_2 * F3_3 - F1_1 * F2_3 * F3_2 - F1_2 * F2_1 * F3_3 + F1_2 * F2_3 * F3_1 + F1_3 * F2_1 * F3_2 - F1_3 * F2_2 * F3_1, 5.0 / 3.0) * (8.0 / 3.0)) + D * pow(F1_1 * F2_2 - F1_2 * F2_1, 2.0) * 2.0;
    */
}