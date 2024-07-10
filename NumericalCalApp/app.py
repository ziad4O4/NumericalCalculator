from tkinter import messagebox
import sympy as sp
import tkinter as tk
from PIL import ImageTk, Image
import os 
import customtkinter as ctk
from tkinter import ttk


def open_window():
    def bisection(func, xl, xu, eps, max_iter=100):
        x = sp.symbols('x')
        f = sp.sympify(func)
        
        if f.subs(x, xl) * f.subs(x, xu) < 0:
            iteration = 0
            xr = 0
            xrold = 0
            error = 100
            iter_result = []
            while error > eps and iteration < max_iter:
                xrold = xr
                xr = (xl + xu) / 2
                error = abs((xr - xrold) / xr) * 100
                if iteration == 0:
                    iter_result.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu), xr,f.subs(x,xr),"____"))
                    print("iteration= ", iteration, "|xl= ", xl, "|f(xl)= ", f.subs(x, xl), "|xu= ", xu, "|f(xu)= ",
                          f.subs(x, xu), "|xr= ", xr, "|f(xr)= ", f.subs(x, xr), "ERROR= ", "")
                else:
                    iter_result.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu), xr,f.subs(x,xr),error))
                if f.subs(x, xl) * f.subs(x, xr) < 0:
                    xu = xr
                else:
                    xl = xr
                iteration += 1
            return xr ,iter_result

    def calculate():
        func = entry_func.get()
        xl   = float(entry_xl.get())
        xu   = float(entry_xu.get())
        eps  = float(entry_eps.get())

        root, result = bisection(func, xl, xu, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f'The Root ="{root}"',font=("bold",16))
 
    window = ctk.CTkToplevel(root)
    window.title("Bisection")
    window.geometry("730x720")
    window.resizable(False,False)
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()

    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 12),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xl = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xl",font=("Arial", 12,), width=170,height=40)
    entry_xl.place(x=40, y=150)

    entry_xu = ctk.CTkEntry(window,placeholder_text="Enter the upper value Xu", font=("Arial", 12), width=170,height=40)
    entry_xu.place(x=245, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 12), width=170, height=40) 
    entry_eps.place(x=450, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" ,command=calculate)
    calculate_button.place(x=50, y=220)
    
    tree = ttk.Treeview(window, columns=("Iteration", "Xl", "F(Xl)", "Xu", "F(Xu)", "Xr", "F(Xr)", 'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=114)
    tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='', text_color='black')
    root_label.place(x=280, y=550)
    
    sol_label = ctk.CTkLabel(window,text="Solution", text_color="black")
    sol_label.place(x=350,y=330)

    window.mainloop()

def false_window():
    def false_position(func, xl, xu, eps, max_iter=100):
        x = sp.symbols('x')
        f = sp.sympify(func)

        if f.subs(x, xl) * f.subs(x, xu) < 0:
            iteration = 0
            xr = 0
            xrold = 0
            error = 100
            iter_result=[]
            while error > eps and iteration < max_iter:
                xrold = xr
                sub = xl - xu
                xr = xu - ((f.subs(x, xu) * sub) / (f.subs(x, xl) - f.subs(x, xu)))
                error = abs((xr - xrold) / xr) * 100
                if iteration == 0:
                    iter_result.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu), xr,f.subs(x,xr),"_____"))
                else:
                    iter_result.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu), xr,f.subs(x,xr),error))
                if f.subs(x, xl) * f.subs(x, xr) < 0:
                    xu = xr
                else:
                    xl = xr
                iteration += 1
            return xr , iter_result

    def calculate():
        func = entry_func.get()
        xl  = float(entry_xl.get())
        xu  = float(entry_xu.get())
        eps = float(entry_eps.get())

        root , result = false_position(func, xl, xu, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f'The Root =" {root} "')

     
    window = ctk.CTkToplevel(root)
    window.title("False Position")
    window.geometry("730x720")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 12),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xl = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xl",font=("Arial", 12), width=170,height=40)
    entry_xl.place(x=40, y=150)

    entry_xu = ctk.CTkEntry(window,placeholder_text="Enter the upper value Xu", font=("Arial", 12), width=170,height=40)
    entry_xu.place(x=245, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 12), width=170, height=40) 
    entry_eps.place(x=450, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
    calculate_button.place(x=50, y=220)
    

    tree = ttk.Treeview(window, columns=("Iteration", "Xl", "F(Xl)", "Xu", "F(Xu)", "Xr", "F(Xr)", 'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=114)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=50, y=550)
    
    sol_label = ctk.CTkLabel(window,text="Solution", text_color="black")
    sol_label.place(x=350,y=330)

    window.mainloop()

def simple_window():
    def simple_fixed(fn, xi, eps, max_iter=100):
        x = sp.symbols("x")
        f = sp.sympify(fn)
        right = 0
        left = 0
        # Get the degrees of the terms in the function
        degrees = [sp.degree(term, gen=x) for term in f.as_ordered_terms()]
        highest_power = max(degrees)
        for term, degree in zip(f.as_ordered_terms(), degrees):
            # Get the term of the highest power at different side of the function
            if degree == highest_power:
                left += (term)
                # Get the coefficients of the terms and change the sign of it in the function
                expr_x = sp.collect(left, x)
                coefficients = sp.Poly(expr_x, x).coeffs()
                left = coefficients[0] * (-1)
            else:
                right += (term)
        # Continue simplifying the function
        g = right / left
        gx = g**(1 / highest_power)
        error = 100.00
        xipls = 0
        iteration = 0
        iter_result = []
        while error > eps and iteration < max_iter:
            if iteration == 0:
                iter_result.append((iteration, xi, f.subs(x,xi),"_____"))
            else:
                xipls = gx.subs(x, xi)
                error = abs((xipls - xi) / xipls) * 100
                xi = xipls
                iter_result.append((iteration, xi, f.subs(x,xi), error))
               
            iteration += 1
        return xi , iter_result
    
    def calculate():
        func = entry_func.get()
        xi   = float(entry_xi.get())
        eps  = float(entry_eps.get())
        root , result = simple_fixed(func, xi, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f'The Root=" {root} "')

    window = ctk.CTkToplevel(root)
    window.title("Simple Fixed Point")
    window.geometry("730x720")
    window.grab_set()

    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 12),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xi = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xi",font=("Arial", 12), width=170,height=40)
    entry_xi.place(x=40, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 12), width=170, height=40) 
    entry_eps.place(x=245, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
    calculate_button.place(x=50, y=220)
    
    tree = ttk.Treeview(window, columns=("Iteration", "Xi", "F(Xi)",'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=220)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=50, y=550)
    
    sol_label = ctk.CTkLabel(window,text="Solution", text_color="black")
    sol_label.place(x=350,y=330)
    
    window.mainloop()

def newton_window():
    def newton(func, xi, eps, max_iter=100):
        x = sp.symbols('x')
        f = sp.sympify(func)
        fdash = sp.diff(f)
        iteration = 0
        xipls = 0
        error = 100.00
        iter_result = []
        iter_result.append((iteration, xi, f.subs(x,xi), fdash.subs(x, xi),"_____"))
        print('Iteration: ', iteration, 'xi: ', xi, 'F(xi): ', f.subs(x, xi), "f'(xi)", fdash.subs(x, xi), 'Error: ', "____")
        iteration += 1
        while error >= eps and iteration < max_iter:
            xipls = xi - (f.subs(x, xi) / fdash.subs(x, xi))
            error = abs((xipls - xi) / xipls) * 100
            xi = xipls
            iter_result.append((iteration, xi, f.subs(x,xi), fdash.subs(x, xi), error))
            print('Iteration: ', iteration, 'xi: ', xi, 'F(xi): ', f.subs(x, xi), "f'(xi)", fdash.subs(x, xi), error)
            iteration += 1
        return xi , iter_result

    def calculate():
        func = entry_func.get()
        xi = float(entry_xi.get())
        eps = float(entry_eps.get())
        root = newton(func, xi, eps, max_iter=100)
        print(f'The ROOT= " {root} "')

        root , result = newton(func, xi, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f'The Root =" {root} "')

    window = ctk.CTkToplevel(root)
    window.title("Newton")
    window.geometry("720x730")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()

    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 12),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xi = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xi",font=("Arial", 12), width=170,height=40)
    entry_xi.place(x=40, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 12), width=170, height=40) 
    entry_eps.place(x=245, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
    calculate_button.place(x=50, y=220)
    
    tree = ttk.Treeview(window, columns=("Iteration", "Xi", "F(Xl)"," F'(Xi) ",'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=175)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=500, y=650)

    sol_label = ctk.CTkLabel(window,text="Solution", text_color="black")
    sol_label.place(x=350,y=330)

    window.mainloop()

def secant_window():
    def secant(func, xmin, xi, eps=1, max_iter=100):
        x = sp.symbols('x')
        f = sp.sympify(func)
        iteration = 0
        xipls = 0
        error = 100
        iter_result = []
        while error > eps and iteration < max_iter:
            xipls = xi - ((f.subs(x, xi) * (xmin - xi)) / (f.subs(x, xmin) - f.subs(x, xi)))
            error = abs((xi - xmin) / xi) * 100
            iter_result.append((iteration, xmin, f.subs(x,xmin),xi,f.subs(x, xi), error))
            print("iteration= ", iteration, "|xi-1= ", xmin, "|f(xi-1)= ", f.subs(x, xmin), "|xi= ", xi, "|f(xi)= ", f.subs(x, xi), "ERROR= ", error, "%")
            iteration += 1
            xmin = xi
            xi = xipls
        return xi , iter_result

    def calculate():
        func = entry_func.get()
        xmin = float(entry_xmin.get())
        xi   = float(entry_xi.get())
        eps  = float(entry_eps.get())
        root = secant(func, xi, eps, max_iter=100)
        print(f'The ROOT =" {root} "')

        root , result = secant(func, xmin,xi, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f"Root= {root}")

     
    window = ctk.CTkToplevel(root)
    window.title("Secant")
    window.geometry("720x730")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 12),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xmin = ctk.CTkEntry(window,placeholder_text="Enter Xmin value ",font=("Arial", 12), width=170,height=40)
    entry_xmin.place(x=40, y=150)

    entry_xi = ctk.CTkEntry(window,placeholder_text="Enter the Xi value ",font=("Arial", 12), width=170,height=40)
    entry_xi.place(x=245, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 12), width=170, height=40) 
    entry_eps.place(x=450, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
    calculate_button.place(x=50, y=220)
    
    tree = ttk.Treeview(window, columns=("Iteration", "Xmin", "F(Xmin)","Xi","F(Xi)",'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=175)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=500, y=650)

    window.mainloop()

def gje_window():
    def gje(matrix):
        def display(matrix):
            for row in matrix:
                print("[", end=" ")
                for elem in row:
                    print("{:5}".format(elem), end=" ")

        m21 = matrix[1][0] / matrix[0][0]
        m31 = matrix[2][0] / matrix[0][0]

        # Rule E2-(m21)E1 = E2
        for j in range(4):
            e2 = matrix[1][j]
            e1 = m21 * matrix[0][j]
            matrix[1][j] = e2 - e1
        new_label.configure(text=f"Step 1")  
        print("step one")
        matrix_label_1_text = "            "
        for row in matrix:
            matrix_label_1_text += "              ".join(map(str, row)) + "\n"
        matrix_label_1.configure(text=f"{matrix_label_1_text}")
        display(matrix)

        # Rule E3-(m31)E1 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e1 = m31 * matrix[0][j]
            matrix[2][j] = e3 - e1
        print("step two")
        matrix_label_2_text = "           "
        for row in matrix:
            matrix_label_2_text += "             ".join(map(str, row)) + "\n"
        matrix_label_2.configure(text=f"{matrix_label_2_text}")
        display(matrix)

        m32 = matrix[2][1] / matrix[1][1]

        # Rule E3-(m32)E2 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e2 = m32 * matrix[1][j]
            matrix[2][j] = e3 - e2
        print("step three")
        matrix_label_3_text = "  "
        for row in matrix:
            matrix_label_3_text += "  ".join(map(str, row)) + "\n"
        matrix_label_3.configure(text=f"{matrix_label_3_text}")
        display(matrix)

        x3 = matrix [2][3] / matrix[2][2]
        x2 = (matrix[1][3] - (matrix[1][2] * x3)) / matrix[1][1]
        x1 = (matrix[0][3] - ((matrix[0][1] * x2) + (matrix[0][2] * x3))) / matrix[0][0]
        
        x3_label.configure(text=f"X3= {x3}")
        x2_label.configure(text=f"X2= {x2}")
        x1_label.configure(text=f"X1= {x1}")
        
    # Main code
        matrix = []
        print("Enter matrix:")
        for i in range(3):
            row = []
            for j in range(4):
                print(f"row{i}{j}: ", end="")
                row.append(float(input()))
            matrix.append(row)

    def calculate():
        matrix = []
        for i in range(3):
            row = []
            for j in range(4):
                value = float(entry_matrix[i][j].get())  # Get the input value and remove leading/trailing whitespace
                row.append(value)
            matrix.append(row)
        gje(matrix)

        # Display the input matrix
        result_text = " "
        for row in matrix:
            result_text += "   ".join(map(str, row)) + "\n"
        display(matrix)
        gje(matrix)  
    
    window = ctk.CTkToplevel(root)
    window.title("Gauss Jordan Elimination")
    window.geometry("720x730")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()

    # Create entry matrix
    entry_matrix = []
    for i in range(3):
        row = []
        for j in range(4):
            entry = ctk.CTkEntry(window, width=45, height=30)
            entry.place(x=100*j+10,y=30*i)
            row.append(entry)
        entry_matrix.append(row)
    for i in range(3):
            for j in range(3):
                xlbl=ctk.CTkLabel(window,text=f"x{j+1}")
                xlbl.place(x=100*j+70,y=30*i)
    
    new_label=ctk.CTkLabel(window,text="")
    new_label.place(x=120,y=220)
    matrix_label_1 = ctk.CTkLabel(window, text="")
    matrix_label_1.place(x=120,y=220)

    matrix_label_2 = ctk.CTkLabel(window, text="")
    matrix_label_2.place(x=350,y=220)
    
    matrix_label_3 = ctk.CTkLabel(window, text="")
    matrix_label_3.place(x=550,y=220)

    x1_label=ctk.CTkLabel(window,text="X1")
    x1_label.place(x=50,y=350)
    x2_label=ctk.CTkLabel(window,text="X2")
    x2_label.place(x=50,y=380)
    x3_label=ctk.CTkLabel(window,text="X3")
    x3_label.place(x=50,y=410)
    

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
    calculate_button.place(x=20, y=140)

    window.mainloop()

def display(matrix):
    for row in matrix:
        print("[", end=" ")
        for elem in row:
            print("{:5}".format(elem), end=" ")
    
def lu_window():
    def lu_decomposition(matrix):
        b = [row[3] for row in matrix]
        print ("b",b ,"\n")
        bvar = "\n".join(str(x) for x in b)
        b_lbl.configure(text=f"B:\n {bvar}")

        #FOR ONE DIMENSIONAL MATRIX
        m21 = matrix[1][0] / matrix[0][0]
        m31 = matrix[2][0] / matrix[0][0]

        # Rule E2-(m21)E1 = E2
        for j in range(4):
            e2 = matrix[1][j]
            e1 = m21 * matrix[0][j]
            matrix[1][j] = e2 - e1

        # Rule E3-(m31)E1 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e1 = m31 * matrix[0][j]
            matrix[2][j] = e3 - e1

        m32 = matrix[2][1] / matrix[1][1]

        # Rule E3-(m32)E2 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e2 = m32 * matrix[1][j]
            matrix[2][j] = e3 - e2
            
            l = [
            [1, 0, 0],
            [m21, 1, 0],
            [m31, m32, 1]
            ]

        l_var = "\n".join(str(x) for x in b )
        l_label.configure(text=f"\n {l_var}")

        u = [[matrix[i][j] for j in range(3)] for i in range(3)]

        # Solve Lc = b
        c = [0] * 3
        c[0] = b[0] / l[0][0]
        c[1] = (b[1] - l[1][0] * c[0]) / l[1][1]
        c[2] = (b[2] - (l[2][0] * c[0] + l[2][1] * c[1])) / l[2][2]
  
        c0_label.configure(text=f"C1= {c[0]}")
        c1_label.configure(text=f"C2= {c[1]}")
        c2_label.configure(text=f"C3= {c[2]}")

        # Solve Ux = c
        x3 = c[2] / u[2][2]
        x2 = (c[1] - u[1][2] * x3) / u[1][1]
        x1 = (c[0] - (u[0][1] * x2 + u[0][2] * x3)) / u[0][0]

        x1_label.configure(text=f"X1= {x1}")
        x2_label.configure(text=f"X2= {x2}")
        x3_label.configure(text=f"X3= {x3}")

    # Main code
        matrix = []
        print("Enter matrix:")
        for i in range(3):
            row = []
            for j in range(4):
                print(f"row{i}{j}: ", end="")
                row.append(float(input()))
            matrix.append(row)
        print("\n")    
        print("the input matrix")    
        display(matrix)    
        lu_decomposition(matrix)  

    def calculate():
        matrix = []
        for i in range(3):
            row = []
            for j in range(4):
                value = float(entry_matrix[i][j].get())  # Get the input value and remove leading/trailing whitespace
                row.append(value)
            matrix.append(row)
        

        # Display the input matrix
        result_text = ""
        for row in matrix:
            result_text += "        ".join(map(str, row)) + "\n"
        matrix_label.configure(text=f"{result_text}")
        
        lu_decomposition(matrix)

    window = ctk.CTkToplevel(root)
    window.title("Lu Decomposition")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()

    # Create entry matrix
    entry_matrix = []
    for i in range(3):
        row = []
        for j in range(4):
            xlbl=ctk.CTkLabel(window,text=f"x{j+1}")
            xlbl.place(x=100*j+70,y=30*i)
            entry = ctk.CTkEntry(window, width=45, height=30)
            entry.place(x=100*j+10,y=30*i)
            row.append(entry)
        entry_matrix.append(row)
     
    matrix_label = ctk.CTkLabel(window, text="", width=100, height=30)
    matrix_label.place(x=200,y=100)

    b_lbl= ctk.CTkLabel(window,text='')
    b_lbl.place(x=40,y=150)

    l_label=ctk.CTkLabel(window,text='')
    l_label.place(x=100 ,y=150)

    c0_label = ctk.CTkLabel(window,text='')
    c0_label.place(x=300,y=150)
    c1_label = ctk.CTkLabel(window,text='')
    c1_label.place(x=300,y=180)
    c2_label = ctk.CTkLabel(window,text='')
    c2_label.place(x=300,y=210)

    x1_label = ctk.CTkLabel(window,text='')
    x1_label.place(x=320,y=150)
    x2_label = ctk.CTkLabel(window,text='')
    x2_label.place(x=310,y=180)
    x3_label = ctk.CTkLabel(window,text='')
    x3_label.place(x=330,y=210)

    calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
    calculate_button.place(x=20, y=140)
########################################################################
def cramer_window():
    def display(matrix):
        for row in matrix:
            print("[", end=" ")
            for elem in row:
                print("{:5}".format(elem), end=" ")
 
         # Copy a matrix
    def copy_matrix(matrix, copy):
        for i in range(3):
            for j in range(3):
                copy[i][j] = matrix[i][j]

    copy = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    def cramer_rule(matrix):
        deta = [0] * 3
        x = [0] * 3
        b = [row[3] for row in matrix]

        # Calculate the delta A
        r0 = matrix[0][0] * ((matrix[1][1] * matrix[2][2]) - (matrix[1][2] * matrix[2][1]))
        r1 = matrix[0][1] * ((matrix[1][0] * matrix[2][2]) - (matrix[1][2] * matrix[2][0]))
        r2 = matrix[0][2] * ((matrix[1][0] * matrix[2][1]) - (matrix[1][1] * matrix[2][0]))
        detAo = r0 - r1 + r2

        detA0_lbl.configure(text=f"Delta A = {detAo} \n")
        # Calculate the delta
        for i in range(3):
            copy_matrix(matrix, copy)
            copy[0][i] = b[0]
            copy[1][i] = b[1]
            copy[2][i] = b[2]

            r0 = copy[0][0] * ((copy[1][1] * copy[2][2]) - (copy[1][2] * copy[2][1]))
            r1 = copy[0][1] * ((copy[1][0] * copy[2][2]) - (copy[1][2] * copy[2][0]))
            r2 = copy[0][2] * ((copy[1][0] * copy[2][1]) - (copy[1][1] * copy[2][0]))
            deta[i] = r0 - r1 + r2

            detAi_lbl.configure(text=detAi_lbl.cget("text") + f"A[{i + 1}] = {deta[i]}\n \n")
        # Calculate the solutions
        for i in range(3):
            x[i] = deta[i] / detAo
            xi_lbl.configure(text=xi_lbl.cget("text") + f" x{i + 1} = {x[i]} ")

    def calculate():
        matrix = []
        for i in range(3):
            row = []
            for j in range(4):
                value = float(entry_matrix[i][j].get())
                row.append(value)
            matrix.append(row)

        result_text = ""
        for row in matrix:
            result_text += "        ".join(map(str, row)) + "\n"
        matrix_lbl.configure(text=result_text)

        cramer_rule(matrix)


    window = ctk.CTk()
    window.title("Cramer")
    window.geometry('730x720')

    # Create entry matrix
    entry_matrix = []
    for i in range(3):
        row = []
        for j in range(4):
            xlbl = ctk.CTkLabel(window, text=f" X {j + 1}")
            xlbl.place(x=100 * j + 70, y=30 * i)
            entry = ctk.CTkEntry(window, width=40, height=1)  
            entry.place(x=100 * j + 10, y=30 * i)
            row.append(entry)
        entry_matrix.append(row)

    calculate_button = ctk.CTkButton(window, text="Calculate", width=220,fg_color="#F2789F", hover_color="#F999B7" , corner_radius=10,command=calculate)
    calculate_button.place(x=100, y=100)

    matrix_lbl = ctk.CTkLabel(window, text="", justify="left")
    matrix_lbl.place(x=150, y=150)

    detA0_lbl = ctk.CTkLabel(window, text='')
    detA0_lbl.place(x=182, y=220)

    detAi_lbl = ctk.CTkLabel(window, text='')
    detAi_lbl.place(x=192, y=270)

    xi_lbl = ctk.CTkLabel(window, text='')
    xi_lbl.place(x=90, y=370)

    root.mainloop()
#####################################################
def partial_window():   
    def rearrange_matrix(matrix):
        if matrix[0][0] < matrix[1][0]:
            temp = matrix[0][0]
            matrix[0][0] = matrix[1][0]
            matrix[1][0] = temp
        if matrix[0][0] < matrix[2][0]:
            temp = matrix[0][0]
            matrix[0][0] = matrix[2][0]
            matrix[2][0] = temp    
        if matrix[1][1] < matrix[2][1]:
            temp = matrix[1][1]
            matrix[1][1] = matrix[2][1]
            matrix[2][1] = temp
        return matrix
    def lu_decomposition_partial(matrix):
        
        b = [row[3] for row in matrix]
        
        m21 = matrix[1][0] / matrix[0][0]
        m31 = matrix[2][0] / matrix[0][0]

        # Rule E2-(m21)E1 = E2
        for j in range(4):
            e2 = matrix[1][j]
            e1 = m21 * matrix[0][j]
            matrix[1][j] = e2 - e1

        # Rule E3-(m31)E1 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e1 = m31 * matrix[0][j]
            matrix[2][j] = e3 - e1

        m32 = matrix[2][1] / matrix[1][1]

        # Rule E3-(m32)E2 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e2 = m32 * matrix[1][j]
            matrix[2][j] = e3 - e2
            
            l = [
            [1, 0, 0],
            [m21, 1, 0],
            [m31, m32, 1]
            ]

        u = [[matrix[i][j] for j in range(3)] for i in range(3)]

        # Solve Lc = b
        c = [0] * 3
        c[0] = b[0] / l[0][0]
        c[1] = (b[1] - l[1][0] * c[0]) / l[1][1]
        c[2] = (b[2] - (l[2][0] * c[0] + l[2][1] * c[1])) / l[2][2]

        c0_label.configure(text=f"C1= {c[0]}")
        c1_label.configure(text=f"C2= {c[1]}")
        c2_label.configure(text=f"C3= {c[2]}")
        # Solve Ux = c
        x  = [0] * 3
        x3 = c[2] / u[2][2]
        x2 = (c[1] - u[1][2] * x[2]) / u[1][1]
        x1 = (c[0] - (u[0][1] * x[1] + u[0][2] * x[2])) / u[0][0]

        x1_label.configure(text=f"X1= {x1}")
        x2_label.configure(text=f"X2= {x2}")
        x3_label.configure(text=f"X3= {x3}")


    def calculate():
        matrix = []
        for i in range(3):
            row = []
            for j in range(4):
                value = float(entry_matrix[i][j].get())  # Get the input value and remove leading/trailing whitespace
                row.append(value)
            matrix.append(row)

        # Display the input matrix
        result_text = ""
        for row in matrix:
            result_text += "        ".join(map(str, row)) + "\n"
        matrix_label.configure(text=f"{result_text}")
        lu_decomposition_partial(matrix)
        

    window = ctk.CTkToplevel(root)
    window.title("Partial")
    window.geometry("700x700")
    window.grab_set()

    # Create entry matrix
    entry_matrix = []
    for i in range(3):
        row = []
        for j in range(4):
            
            entry = ctk.CTkEntry(window, width=45, height=30)
            entry.place(x=100*j+10,y=30*i)
            row.append(entry)
        entry_matrix.append(row)

    for i in range(3):
        for j in range(3):
            xlbl=ctk.CTkLabel(window,text=f"x{j+1}")
            xlbl.place(x=100*j+70,y=30*i)
 
        matrix_label = ctk.CTkLabel(window, text="", width=100, height=30)
        matrix_label.place(x=2,y=100)

        x1_label = ctk.CTkLabel(window,text='')
        x1_label.place(x=400,y=200)
        x2_label = ctk.CTkLabel(window,text='')
        x2_label.place(x=400,y=200)
        x3_label = ctk.CTkLabel(window,text='')
        x3_label.place(x=400,y=250)
        
        c0_label = ctk.CTkLabel(window,text='')
        c0_label.place(x=300,y=150)
        c1_label = ctk.CTkLabel(window,text='')
        c1_label.place(x=300,y=180)
        c2_label = ctk.CTkLabel(window,text='')
        c2_label.place(x=300,y=210)
        
        calculate_button = ctk.CTkButton(window, text="Calculate",fg_color="#F2789F", hover_color="#F999B7" , command=calculate)
        calculate_button.place(x=20, y=140)

    window.mainloop()

def golden_window():
    def determinMax(func,xl,xu,stop_iter):
        x = sp.symbols('x')
        f = sp.sympify(func)
        gr = ((5**0.5)- 1) / 2  # Golden ratio
        d = gr * (xu - xl)  # Initial distance
        x1 = xu - d  # First trial point
        x2 = xl + d  # Second trial point
        iteration = 1
        iter_res=[]
        iter_res.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu),x1,f.subs(x,x1),x2,f.subs(x,x2),d))
        while iteration < stop_iter :
            iteration +=1
            if(f.subs(x,x1)<f.subs(x,x2)) :
                xu=x1
                x1=x2
            else :
                xl=x2
                x2=x1
            d = gr * (xu - xl)  # Initial distance
            x1 = xu - d  # First trial point
            x2 = xl + d  # Second trial point
            
            iter_res.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu),x1,f.subs(x,x1),x2,f.subs(x,x2),d))
        if (x1 > x2):
            max=x1
        else :
            max=x2    
        return iter_res, max

    def determinMin(func,xl,xu,site):
        x = sp.symbols('x')
        f = sp.sympify(func)
        gr = ((5**0.5)- 1) / 2  # Golden ratio
        d = gr * (xu - xl)  # Initial distance
        x1 = xu - d  # First trial point
        x2 = xl + d  # Second trial point
        iteration = 1
        iter_res=[]
        iter_res.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu),x1,f.subs(x,x1),x2,f.subs(x,x2),d))
        while iteration < site :
            iteration +=1
            if(f.subs(x,x1)<f.subs(x,x2)) :
                xl=x2
                x2=x1
                
            else :
                xu=x1
                x1=x2
            d = gr * (xu - xl)  # Initial distance
            x1 = xu - d  # First trial point
            x2 = xl + d  # Second trial point
            
            iter_res.append((iteration,xl,f.subs(x,xl),xu,f.subs(x,xu),x1,f.subs(x,x1),x2,f.subs(x,x2),d))
        if (x1 < x2):
            min=x1
        else :
            min=x2    
        return iter_res, min

    #MAIN CODE
    # Input function from the user
    def max_calculate_root():
        fn = function_entry.get()
        xl = float(xl_entry.get())
        xu = float(xu_entry.get())
        stop_iter = float(stop_entry.get())

        iteration_results, root = determinMax(fn, xl, xu, stop_iter)
        for result in iteration_results:
            tree.insert("", "end", values=result)
        root_label.configure(text=f"ROOT: {root} %")


    def min_calculate_root():
        fn = function_entry.get()
        xl = float(xl_entry.get())
        xu = float(xu_entry.get())
        stop_iter = float(stop_entry.get())

        iteration_results, root = determinMin(fn, xl, xu, stop_iter)
        for result in iteration_results:
            tree.insert("", "end", values=result)
        root_label.configure(text=f"ROOT: {root} %")

    window = ctk.CTk()
    window.geometry('730×720')
    window.title("Golden Method")
    window.resizable(False, False)

    function_entry = ctk.CTkEntry(window, placeholder_text="Enter the function",width=200)
    function_entry.place(x=180, y=10)

    # XL entry
    xl_entry = ctk.CTkEntry(window, placeholder_text='Enter the value of Xl',width=200)
    xl_entry.place(x=180, y=50)

    # XU entry
    xu_entry = ctk.CTkEntry(window, placeholder_text='Enter the value of Xu',width=200)
    xu_entry.place(x=180, y=90)

    # Error entry
    stop_entry = ctk.CTkEntry(window, placeholder_text='stop iteration',width=200)
    stop_entry.place(x=180, y=130)

    # Buttons
    max_calculate_button = ctk.CTkButton(window, text="Calculate MAX", width=250, command=max_calculate_root)
    max_calculate_button.place(x=840, y=20)

    min_calculate_button = ctk.CTkButton(window, text="Calculate MIN", width=250, command=min_calculate_root)
    min_calculate_button.place(x=840, y=70)

    tree = ttk.Treeview(window, columns=('Iteration', 'Xl', 'f(XL)', 'Xu', 'f(Xu)', 'X1', 'f(X1)', 'X2', 'f(X2)', 'd'), show='headings')
    tree.column("#0", width=30)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=104)
    tree.place(x=10, y=200)

    # root
    root_label = ctk.CTkLabel(window, text='', font=("bold", 19))
    root_label.place(x=10,y=500)

    window.mainloop()

root = ctk.CTk()
root.geometry("800x700")
root.title("Numerical Calculator")

ctk.set_default_color_theme("green")
ctk.set_appearance_mode("light")    

# Load the background image
img_path = os.path.join(os.path.dirname(__file__), r'D:\\code\\python\\Numerical_Cal_App\\phhhh.jpeg')
img = ctk.CTkImage(light_image=Image.open(img_path), size=(1050, 5700))

# Create a label to display the background image
img_lbl = ctk.CTkLabel(root, image=img, text='')
img_lbl.place(x=0, y=0, relwidth=1, relheight=1)

method_label = ctk.CTkLabel(root, text="Choose the method you want ",font=("Arial", 24, "bold"), anchor = ctk.CENTER)
method_label.pack(pady=50)
########
bisection_button = ctk.CTkButton(root,text="Bisection", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER, command=open_window)
bisection_button.pack(pady=10)
########
false_position_button = ctk.CTkButton(root, text="False Position", height=30, width=240,fg_color="#F999B7", hover_color="#F9C5D5",corner_radius=8, anchor=ctk.CENTER,
                                          command=false_window)
false_position_button.pack(pady=10)
#########
simple_fixed_button = ctk.CTkButton(root, text="Simple Fixed", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5",corner_radius=8, anchor=ctk.CENTER,
                                        command=simple_window)
simple_fixed_button.pack(pady=10)
##########
newton_button = ctk.CTkButton(root, text="Newton", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5",corner_radius=8, anchor=ctk.CENTER,
                                  command=newton_window)
newton_button.pack(pady=10)
########
secant_button = ctk.CTkButton(root, text="Secant", height=30, width=240,fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER,
                                  command=secant_window)
secant_button.pack(pady=10)
#########
guass_jordan_button = ctk.CTkButton(root, text="Guass Jordan Elimination", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER,
                                        command=gje_window)
guass_jordan_button.pack(pady=10)
##########
lu_button = ctk.CTkButton(root, text="Lu Decomposition", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER,
                                        command=lu_window)
lu_button.pack(pady=10)
##########
cramer_button = ctk.CTkButton(root, text="Cramer", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER,
                                        command=cramer_window)
cramer_button.pack(pady=10)
##########
partial_button = ctk.CTkButton(root, text="Partial", height=30, width=240, fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER, command=partial_window)
partial_button.pack(pady=10)
#########
golden_button = ctk.CTkButton(root, text="Golden Section", height=30, width=240,  fg_color="#F999B7", hover_color="#F9C5D5", corner_radius=8, anchor=ctk.CENTER,
                                       command=golden_window)
golden_button.pack(pady=10)

root.mainloop()

