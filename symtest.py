

from sympy import *;
from sympy import symbols
from sympy.physics.vector import *;
from sympy.vector import matrix_to_vector

C = ReferenceFrame('C')
t,s = symbols('t s');
pprint_use_unicode(false);

def speed(r_t,a,b):
    v_t = r_t.diff(t, C)

    print("Velocity \n");
    pprint(v_t ,use_unicode=false);

    print("\n Speed \n")
    speed = v_t.magnitude().simplify()
    pprint(speed,use_unicode = false);

    print("\n Length of trajectory \n")
    pprint(integrate(speed,(t,a,b)))

def torsion(r_t):
    v = r_t.diff(t,C)
    a = v.diff(t,C)
    aprime = a.diff(t,C)
    unit_tangent = v/(v.magnitude());

    normal = unit_tangent.diff(t,C)
    unit_normal = normal / normal.magnitude();

    print("\n Binormal \n")
    pprint(cross(unit_tangent,unit_normal).simplify())

    numerator = (cross(v,a)).dot(aprime)

    denominator = ((cross(v,a)).magnitude())**2

    tor = (numerator/denominator).simplify()

    print("\n Torsion \n")
    pprint(tor)

def curvature(r_t):
    v = r_t.diff(t,C);
    a = v.diff(t,C);

    print("\n Curvature \n")

    pprint(((cross(v,a).magnitude())/(v.magnitude())**3).simplify());

def arclength_check(r_t):
    v = r_t.diff(t,C).simplify();
    if v.magnitude().simplify() == 1:
        print("\nThe curve uses arc length as a parameter\n");

    else:
        expr = (integrate(v.magnitude(),t).simplify())-s;
        expr = solve(expr,t)
        r_s = r_t.subs(t,expr[0]);

        print("the reparameterized curve")
        pprint(r_s);

def function_builder(components):
    sympyComponents = [0,0,0]
    for i in range(len(components)):
        sympyComponents[i] = parse_expr(components[i]);

    vector = sympyComponents[0] * C.x + sympyComponents[1] * C.y + sympyComponents[2] * C.z;

    return vector;

def choices(choice,r_t,bounds):
    if choice == 1:
        speed(r_t, parse_expr(bounds[0]), parse_expr(bounds[1]));

    elif choice == 2:
        curvature(r_t);

    elif choice == 3:
        torsion(r_t)

    elif choice == 4:
        arclength_check(r_t);

def main():
    while True:

        print("\n1. Speed and Trajectory of Curve")
        print("2. Curvature of Curve");
        print("3. Torsion of Curve");
        print("4. Arclength Reparameterization of Curve\n");
        print("\npress 1 through 4: \n");
        choice = int(input());

        components = [0,0,0];

        for i in range(0,3):
            print("Enter vector component: ")
            user_input = str(input());
            components[i] = user_input;

        bounds = ["0","0"]

        if choice == 1:
            print("Enter bounds as comma separated list: ");
            temp = input();
            bounds = temp.split(",");

        r_t = function_builder(components);

        choices(choice,r_t,bounds)

main()


































