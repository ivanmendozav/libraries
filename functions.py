from sympy import symbols, Eq, Matrix, solve, Rational, latex, Symbol
import numpy as np

def valores_propios(A):
  A = Matrix(A)
  L = Symbol("\lambda")
  n = A.shape[0]

  A2 =L*Matrix.eye(n)-A
  A2.det()

  Px = A2.det()
  print("\nPolinomio caracteristico:\n")

  print("\nValores propios:\n")
  sols= solve(Eq(Px,0),L) 

  for i in range(len(sols)):
    A3 = A2.subs(L,sols[i])
    sol = A3.gauss_jordan_solve(Matrix.zeros(1,n).T)
    print(f"\nVectores propios:\n")

def Eigen2Matrix3x3(valores, vector1, vector2, vector3):
  """
  Args:
    valores:Eigen valores
    vector1:Eigen vector para lambda1
    vector2:Eigen vector para lambda2
    vector3:Eigen vector para lambda3
  Returns:Matriz transformacion

  """
  a11, a12, a13, a21, a22, a23, a31, a32, a33 = symbols('a11 a12 a13 a21 a22 a23 a31 a32 a33')

  l1 = valores[0]
  l2 = valores[1]
  l3 = valores[2]

  # Definir sistema de ecuaciones
  Eq1 = Eq((a11 - l1)*vector1[0] + a12*vector1[1] + a13*vector1[2], 0)
  Eq2 = Eq(a21*vector1[0] + (a22 - l1)*vector1[1] + a23*vector1[2], 0)
  Eq3 = Eq(a31*vector1[0] + a32*vector1[1] + (a33 - l1)*vector1[2], 0)

  Eq4 = Eq((a11 - l2)*vector2[0] + a12*vector2[1] + a13*vector2[2], 0)
  Eq5 = Eq(a21*vector2[0] + (a22 - l2)*vector2[1] + a23*vector2[2], 0)
  Eq6 = Eq(a31*vector2[0] + a32*vector2[1] + (a33 - l2)*vector2[2], 0)

  Eq7 = Eq((a11 - l3)*vector3[0] + a12*vector3[1] + a13*vector3[2], 0)
  Eq8 = Eq(a21*vector3[0] + (a22 - l3)*vector3[1] + a23*vector3[2], 0)
  Eq9 = Eq(a31*vector3[0] + a32*vector3[1] + (a33 - l3)*vector3[2], 0)

  # Resolver para elementos de la matriz
  B = solve((Eq1, Eq2, Eq3, Eq4, Eq5, Eq6, Eq7, Eq8, Eq9))
  A = Matrix([ [B[a11],B[a12],B[a13]],[B[a21],B[a22],B[a23]],[B[a31],B[a32],B[a33]] ])

  # Comprobacion
  valores_propios(A)

  return A

# Prueba
if __name__ == "__main__":
    vector1 = Matrix([2, -1, 1]).T
    vector2 = Matrix([1, 0, -1]).T
    vector3 = Matrix([2, 2, 1]).T
    valores = Matrix([-1, -2, -3]).T
    Eigen2Matrix3x3(valores, vector1, vector2, vector3)