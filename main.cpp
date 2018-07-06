#include <fstream>
#include <iostream>
#include <vector>
#include "simple_matrix.hpp"

using namespace std;


SMatrix<float> add(SMatrix<float> Z1, SMatrix<float> Z2) {
    SMatrix<float> Z3(Z1.size());
  
    int const H1 = (Z1.size());
    if (Z1.size() == 1) {
        Z3(0, 0) = Z1(0, 0) + Z2(0, 0);
        return Z3;
    } else {
        for (int i = 0; i < H1; i++) {
            for (int j = 0; j < H1; j++) {
                Z3(i, j) = Z1(i, j) + Z2(i, j);
            }
        }

        return Z3;
    }

}


SMatrix<float> sub(SMatrix<float> Z1, SMatrix<float> Z2) {
    SMatrix<float> Z3(Z1.size());

    int const H1 = (Z1.size());
    if (Z1.size() == 1) {
        Z3(0, 0) = Z1(0, 0) - Z2(0, 0);
        return Z3;
    } else {
        for (int i = 0; i < H1; i++) {
            for (int j = 0; j < H1; j++) {
                Z3(i, j) = Z1(i, j) - Z2(i, j);
            }
        }

        return Z3;
    }

}

SMatrix<float> recmul(SMatrix<float> A, SMatrix<float> B) {
    SMatrix<float> C(A.size());
    if (A.size() == 1) {
        C(0, 0) = A(0, 0) * B(0, 0);
    } else {
   const int N=A.size();

        const int H = (A.size() / 2);

        SMatrix<float> A11(H), A12(H), A21(H), A22(H), B11(H), B12(H), B21(H), B22(H);

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < H; j++) {
                A11(i, j) = A(i, j);
                B11(i, j) = B(i, j);
            }
        }

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < H; j++) {
                A12(i, j) = A(i, j+H);
                B12(i, j) = B(i, j+H);
            }
        }

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < H; j++) {
                A21(i, j) = A(i+H, j);
                B21(i, j) = B(i+H, j);
            }
        }

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < H; j++) {
                A22(i, j) = A(i + H, j + H);
                B22(i, j) = B(i + H, j + H);
            }
        }

        SMatrix<float> S1(H), S2(H), S3(H), S4(H), S5(H), S6(H), S7(H), S8(H), S9(H), S10(H);
        S1 = sub(B12, B22);
        S2 = add(A11, A12);
        S3 = add(A21, A22);
        S4 = sub(B21, B11);
        S5 = add(A11, A22);
        S6 = add(B11, B22);
        S7 = sub(A12, A22);
        S8 = add(B21, B22);
        S9 = sub(A11, A21);
        S10 = add(B11, B12);
        SMatrix<float> X1(H), X2(H), X3(H), X4(H), X5(H), X6(H), X7(H);
       
             if (N > 2) {
            X1 = recmul(A11, S1);
            X2 = recmul(S2, B22);
            X3 = recmul(S3, B11);
            X4 = recmul(A22, S4);
            X5 = recmul(S5, S6);
            X6 = recmul(S7, S8);
            X7 = recmul(S9, S10);
        } else {
            X1(0, 0) = A11(0, 0) * S1(0, 0);
            X2(0, 0) = S2(0, 0) * B22(0, 0);
            X3(0, 0) = S3(0, 0) * B11(0, 0);
            X4(0, 0) = A22(0, 0) * S4(0, 0);
            X5(0, 0) = S5(0, 0) * S6(0, 0);
            X6(0, 0) = S7(0, 0) * S8(0, 0);
            X7(0, 0) = S9(0, 0) * S10(0, 0);
        }

            X1 = recmul(A11, S1);
            X2 = recmul(S2, B22);
            X3 = recmul(S3, B11);
            X4 = recmul(A22, S4);
            X5 = recmul(S5, S6);
            X6 = recmul(S7, S8);
            X7 = recmul(S9, S10);

        SMatrix<float> C11(H), C12(H), C21(H), C22(H);


      



        C11 = add(sub(add(X5,X4),X2),X6);
        C12 = add(X1,X2);  // function
        C21 = add(X3,X4);
        C22 = sub(sub(add(X5,X1),X3),X7);

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < H; j++) {


                C(i, j) = C11(i, j);
                C(i, j + H) = C12(i, j);
                C(i + H, j) = C21(i, j);
                C(i + H, j + H) = C22(i, j);

            }
        }


    }
    return C;

}

int main(int argc, char *argv[]) {


    ifstream in(argv[1]);

    if (!in) {
        cout << "Cannot open  input file;" << endl;
        return -1;
    }

    ofstream out(argv[2]);

    if (!out) {
        cout << "Cannot open  output file;" << endl;
        return -1;
    }
    int datasets;
    in >> datasets;


    for (int i = 0; i < datasets; i++) {
        int N;
        in >> N;


        SMatrix<float> A(N), B(N), C(N);


        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                in>>A.operator()(j, k);
            }
        }
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                in >> B.operator()(j, k);
            }
        }

        C = recmul(A, B);


        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                out << C(i, j) << " ";
            }
        }

        out << endl;
    }

    return 0;
}





