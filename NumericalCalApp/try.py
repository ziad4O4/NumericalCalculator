from tkinter import messagebox
import numpy as np
import sympy as sp
import tkinter as tk
# from PIL import ImageTk, Image
import customtkinter as ctk
from tkinter import ttk


# from scipy import optimize


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
        root_label.configure(text=f"Root= {root}")

     
    window = ctk.CTkToplevel(root)
    window.title("Bisection")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 10),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xl = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xl",font=("Arial", 10), width=170,height=40)
    entry_xl.place(x=40, y=150)

    entry_xu = ctk.CTkEntry(window,placeholder_text="Enter the upper value Xu", font=("Arial", 10), width=170,height=40)
    entry_xu.place(x=245, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 10), width=170, height=40) 
    entry_eps.place(x=450, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate", command=calculate)
    calculate_button.place(x=290, y=270)
    

    tree = ttk.Treeview(window, columns=("Iteration", "Xl", "F(Xl)", "Xu", "F(Xu)", "Xr", "F(Xr)", 'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=110)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='', text_color='black')
    root_label.place(x=500, y=650)
    
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
        root_label.configure(text=f"Root= {root}")

     
    window = ctk.CTkToplevel(root)
    window.title("False Position")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 10),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xl = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xl",font=("Arial", 10), width=170,height=40)
    entry_xl.place(x=40, y=150)

    entry_xu = ctk.CTkEntry(window,placeholder_text="Enter the upper value Xu", font=("Arial", 10), width=170,height=40)
    entry_xu.place(x=245, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 10), width=170, height=40) 
    entry_eps.place(x=450, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate", command=calculate)
    calculate_button.place(x=290, y=270)
    

    tree = ttk.Treeview(window, columns=("Iteration", "Xl", "F(Xl)", "Xu", "F(Xu)", "Xr", "F(Xr)", 'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=110)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=500, y=650)
    
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
        root_label.configure(text=f"Root= {root}")

     
    window = ctk.CTkToplevel(root)
    window.title("Simple Fixed Point")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 10),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xi = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xi",font=("Arial", 10), width=170,height=40)
    entry_xi.place(x=40, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 10), width=170, height=40) 
    entry_eps.place(x=245, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate", command=calculate)
    calculate_button.place(x=290, y=270)
    

    tree = ttk.Treeview(window, columns=("Iteration", "Xi", "F(Xi)",'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=218)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=500, y=650)

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
        print(f"ROOT= {root}")

        root , result = newton(func, xi, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f"Root= {root}")

     
    window = ctk.CTkToplevel(root)
    window.title("Newton")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 10),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xi = ctk.CTkEntry(window,placeholder_text="Enter the lower value Xi",font=("Arial", 10), width=170,height=40)
    entry_xi.place(x=40, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 10), width=170, height=40) 
    entry_eps.place(x=245, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate", command=calculate)
    calculate_button.place(x=290, y=270)
    

    tree = ttk.Treeview(window, columns=("Iteration", "Xi", "F(Xl)"," F'(Xi) ",'Error (%)'), show='headings')
    tree.column("#0", width=50)
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=175)
        tree.place(x=1, y=450)

    root_label = ctk.CTkLabel(window, text='',text_color='black')
    root_label.place(x=500, y=650)

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
        print(f"ROOT= {root}")

        root , result = secant(func, xmin,xi, eps, max_iter=100)
        for rst in result:
            tree.insert('',tk.END,values=rst)
        root_label.configure(text=f"Root= {root}")

     
    window = ctk.CTkToplevel(root)
    window.title("Secant")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()


    entry_func = ctk.CTkEntry(window,placeholder_text="Enter the function f(x)",font=("Arial", 10),width=580,height=40)
    entry_func.place(x=40, y=100)

    entry_xmin = ctk.CTkEntry(window,placeholder_text="Enter Xmin value ",font=("Arial", 10), width=170,height=40)
    entry_xmin.place(x=40, y=150)

    entry_xi = ctk.CTkEntry(window,placeholder_text="Enter the Xi value ",font=("Arial", 10), width=170,height=40)
    entry_xi.place(x=245, y=150)

    entry_eps = ctk.CTkEntry(window,placeholder_text="Enter the epsilon value (Eps)", font=("Arial", 10), width=170, height=40) 
    entry_eps.place(x=450, y=150)

    calculate_button = ctk.CTkButton(window, text="Calculate", command=calculate)
    calculate_button.place(x=290, y=270)
    

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
                print("]")
            print("\n")

        m21 = matrix[1][0] / matrix[0][0]
        m31 = matrix[2][0] / matrix[0][0]

        # Rule E2-(m21)E1 = E2
        for j in range(4):
            e2 = matrix[1][j]
            e1 = m21 * matrix[0][j]
            matrix[1][j] = e2 - e1
        print("step one")
        display(matrix)

        # Rule E3-(m31)E1 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e1 = m31 * matrix[0][j]
            matrix[2][j] = e3 - e1
        print("step two")
        display(matrix)

        m32 = matrix[2][1] / matrix[1][1]

        # Rule E3-(m32)E2 = E3
        for j in range(4):
            e3 = matrix[2][j]
            e2 = m32 * matrix[1][j]
            matrix[2][j] = e3 - e2
        print("step three")
        display(matrix)

        x3 = matrix[2][3] / matrix[2][2]
        x2 = (matrix[1][3] - (matrix[1][2] * x3)) / matrix[1][1]
        x1 = (matrix[0][3] - ((matrix[0][1] * x2) + (matrix[0][2] * x3))) / matrix[0][0]

        return x1, x2, x3

    def calculate():
        matrix = []
        for i in range(3):
            row = []
            for j in range(4):
                value = float(entry_matrix[i][j].get())  # Get the input value
                row.append(value)
            matrix.append(row)

        x1, x2, x3 = gje(matrix)

        x1label.configure(text=f"X1= {x1}")
        x1label.place(x=100,y=100)
        x2label.configure(text=f"X2= {x2}")
        x2label.place(x=200,y=200)
        x3label.configure(text=f"X3= {x3}")
        x3label.place(x=300,y=300)

        # Display the input matrix
        result_text = ""
        for row in matrix:
            result_text += "        ".join(map(str, row)) + "\n"
        matrix_label.configure(text=f"{result_text}")

    # Create Tkinter application
    root = ctk.CTk()
    root.title("Gauss Jordan Elimination")
    root.geometry("700x700")

    # Create entry matrix
    entry_matrix = []
    for i in range(3):
        row = []
        for j in range(4):
            xlbl = ctk.CTkLabel(window, text=f"x{j+1}")
            xlbl.place(x=100*j+70, y=30*i)
            entry = ctk.CTkEntry(window, width=5, height=1)
            entry.place(x=100*j+10, y=30*i)
            row.append(entry)
        entry_matrix.append(row)
    
    window = ctk.CTkToplevel(root)
    window.title("Guass jordan Elimination")
    window.geometry("700x700")
    ctk.set_colors(window, background="white", foreground="white")
    window.grab_set()
    
    matrix_label = ctk.CTkLabel(window, text="", width=100, height=10)
    matrix_label.place(x=2, y=100)

    x1label = ctk.CTkLabel(window, text="")
    x1label.place(x=10, y=250)
    x2label = ctk.CTkLabel(window, text="")
    x2label.place(x=10, y=270)
    x3label = ctk.CTkLabel(window, text="")
    x3label.place(x=10, y=290)

    calculate_button = ctk.CTkButton(window, text="Calculate", command=calculate)
    calculate_button.place(x=290, y=270)

    window.mainloop()

    
root = ctk.CTk()
root.geometry("800x700")
root.title("Numerical Calculator")

ctk.set_default_color_theme("green")
ctk.set_appearance_mode("light")    

method_label = ctk.CTkLabel(root, text="Choose your method",font=("Helvetica", 24, "bold"), anchor = ctk.CENTER)
method_label.pack(pady=50)
########
bisection_button = ctk.CTkButton(root, text="Bisection", height=20, width=180, hover_color="green", corner_radius=8, anchor=ctk.CENTER, command=open_window)
bisection_button.pack(pady=10)
########
false_position_button = ctk.CTkButton(root, text="False Position", height=20, width=180, hover_color="green",corner_radius=8, anchor=ctk.CENTER,
                                          command=false_window)
false_position_button.pack(pady=10)
#########
simple_fixed_button = ctk.CTkButton(root, text="Simple Fixed", height=20, width=180, hover_color="green", corner_radius=8, anchor=ctk.CENTER,
                                        command=simple_window)
simple_fixed_button.pack(pady=10)
##########
newton_button = ctk.CTkButton(root, text="Newton", height=20, width=180, hover_color="green",corner_radius=8, anchor=ctk.CENTER,
                                  command=newton_window)
newton_button.pack(pady=10)
########
secant_button = ctk.CTkButton(root, text="Secant", height=20, width=180, hover_color="green", corner_radius=8, anchor=ctk.CENTER,
                                  command=secant_window)
secant_button.pack(pady=10)
#########
guass_jordan_button = ctk.CTkButton(root, text="Guass Jordan Elimination", height=20, width=180, hover_color="green", corner_radius=8, anchor=ctk.CENTER,
                                        command=gje_window)
guass_jordan_button.pack(pady=10)


root.mainloop()

              
