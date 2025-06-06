import tkinter as tk
from tkinter import ttk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

root = tk.Tk()
check_var = tk.BooleanVar()

class MathResultsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Результаты вычислений")
        self.math_font = ('Arial', 12)
        
        self.header_frame = ttk.Frame(root, padding="10")
        self.main_frame = ttk.Frame(root, padding="10")
        self.table_frame = ttk.Frame(root)
        
        self.header_frame.pack(fill=tk.X)
        ttk.Separator(root, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        self.create_header()
        self.load_initial_data()
        
    def create_header(self):
        try:
            with open("./output/parameters.txt") as f:
                params = [f.readline().strip() for _ in range(4)]
                
            ttk.Label(self.header_frame, text=f"n = {params[0]}", 
                     font=self.math_font).grid(row=0, column=0, padx=5)
            ttk.Label(self.header_frame, text=f"m = {params[1]}", 
                     font=self.math_font).grid(row=0, column=1, padx=5)
            eps = params[2]
            ttk.Label(self.header_frame, text=f"ε = {eps}", 
                     font=self.math_font).grid(row=0, column=2, padx=5)
            ttk.Label(self.header_frame, text=f"N_max = {params[3]}", 
                     font=self.math_font).grid(row=0, column=3, padx=5)
            
            ttk.Radiobutton(self.header_frame, text="Тест", 
                          variable=check_var, value=False,
                          command=self.update_display).grid(row=0, column=4, padx=5)
            ttk.Radiobutton(self.header_frame, text="Основная", 
                          variable=check_var, value=True,
                          command=self.update_display).grid(row=0, column=5, padx=5)
            
        except Exception as e:
            print(f"Ошибка загрузки параметров: {e}")

    def load_initial_data(self):
        check_var.set(False)
        self.update_display()

    def update_display(self):
        for widget in self.main_frame.winfo_children():
            widget.destroy()
            
        for widget in self.table_frame.winfo_children():
            widget.destroy()
            
        self.table_frame.pack_forget()

        if check_var.get(): 
            self.show_main_results()
        else:  
            self.show_test_results()

    def show_main_results(self):
        try:
            with open("./output/main_results.txt") as f:
                results = [f.readline().strip() for _ in range(11)]
                
            main_info = ttk.Frame(self.main_frame)
            main_info.pack(fill=tk.X, pady=10)
            
            ttk.Label(main_info, text="Основная сетка (n, m)", 
                     font=self.math_font).grid(row=0, column=0, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Итераций: {results[0]}", 
                     font=self.math_font).grid(row=1, column=0, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Норма нулевой невязки: {float(results[7]):.2e}", 
                     font=self.math_font).grid(row=2, column=0, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Норма конечной невязки: {float(results[9]):.2e}", 
                     font=self.math_font).grid(row=3, column=0, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Точность: {float(results[1]):.2e}", 
                     font=self.math_font).grid(row=4, column=0, sticky=tk.W, padx=5)
            
            ttk.Label(main_info, text="Удвоенная сетка (2n, 2m)", 
                     font=self.math_font).grid(row=0, column=1, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Итераций: {results[2]}", 
                     font=self.math_font).grid(row=1, column=1, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Норма нулевой невязки: {float(results[8]):.2e}", 
                     font=self.math_font).grid(row=2, column=1, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Норма конечной невязки: {float(results[10]):.2e}", 
                     font=self.math_font).grid(row=3, column=1, sticky=tk.W, padx=5)
            ttk.Label(main_info, text=f"Точность: {float(results[3]):.2e}", 
                     font=self.math_font).grid(row=4, column=1, sticky=tk.W, padx=5)
            
            ttk.Label(main_info, 
                     text=f"max|u-u°|: {float(results[4]):.2e} при x={float(results[5]):.3f}, y={float(results[6]):.3f}", 
                     font=self.math_font).grid(row=5, column=0, columnspan=2, sticky=tk.W, padx=5)
            
            df_v1 = pd.read_csv("./output/v1_main.csv")
            df_v = pd.read_csv("./output/v2_main.csv")
            df_uv = pd.read_csv("./output/uv_main.csv")
            
            self.show_tables(df_v1, df_v, df_uv, 1, "v(n, m)", "v(2n, 2m)", "|v - v2|")
            
        except Exception as e:
            print(f"Ошибка загрузки основных результатов: {e}")

    def show_test_results(self):
        try:
            with open("./output/test_results.txt") as f:
                results = [f.readline().strip() for _ in range(11)]
                
            ttk.Label(self.main_frame, text=f"Итераций: {results[0]}", 
                     font=self.math_font).pack(anchor=tk.W)
            ttk.Label(self.main_frame, text=f"Норма нулевой невязки: {float(results[1][:4] + results[1][12:]):.2e}", 
                     font=self.math_font).pack(anchor=tk.W)
            ttk.Label(self.main_frame, text=f"Точность: {float(results[2])}", 
                     font=self.math_font).pack(anchor=tk.W)
            ttk.Label(self.main_frame, text=f"Норма конечной невязки: {float(results[3]):.2e}", 
                     font=self.math_font).pack(anchor=tk.W)
            ttk.Label(self.main_frame, 
                     text=f"max|u-u°|: {float(results[4]):.2e} при x={float(results[6]):.3f}, y={float(results[5]):.3f}", 
                     font=self.math_font).pack(anchor=tk.W)
            
            df_u = pd.read_csv("./output/u_test.csv")
            df_v = pd.read_csv("./output/v_test.csv")
            df_uv = pd.read_csv("./output/uv_test.csv")
            
            self.show_tables(df_u, df_v, df_uv, 0, "u", "v", "|u-v|")

                     
        except Exception as e:
            print(f"Ошибка загрузки тестовых результатов: {e}")

    def show_tables(self, df_u, df_v, df_uv, fl, n1, n2, n3):
        try:
            self.table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
            
            notebook = ttk.Notebook(self.table_frame)
            notebook.pack(fill=tk.BOTH, expand=True)
            
            with open("./output/parameters.txt") as f:
                params = [int(f.readline().strip()) for _ in range(2)]
                
            h = 2 / params[0]
            k = 1 / params[1]            
            
            self.create_table_tab(notebook, n1, df_u, h, k)
            if (fl == 0):
                self.create_table_tab(notebook, n2, df_v, h, k)
            else:
                self.create_table_tab(notebook, n2, df_v, h/2, k/2)

            self.create_table_tab(notebook, n3, df_uv, h, k)
            
            if (fl == 0):
                self.create_plots_tab(notebook, df_u.values, df_v.values, df_uv.values, h, k)
                print(df_u.values)
            else:
                self.create_plots_tab1(notebook, df_u.values, df_v.values, df_uv.values, h, k)
                print(df_u.values)

        except Exception as e:
            print(f"Ошибка загрузки таблиц: {e}")

    def create_table_tab(self, notebook, title, df, h, k):
        tab = ttk.Frame(notebook)
        notebook.add(tab, text=title)
        
        tree = ttk.Treeview(tab, show="headings")
        
        x_coords = [round(i * h, 4) for i in range(df.shape[1])]
        tree["columns"] = ["y/x"] + x_coords
        
        tree.heading("y/x", text="y \ x")
        tree.column("y/x", width=80, anchor='center')
        
        for col in x_coords:
            tree.heading(col, text=str(col))
            tree.column(col, width=80, anchor='center')
        
        y_coords = [round(j * k, 4) for j in range(df.shape[0])]
        for idx, row in df.iterrows():
            tree.insert("", "end", values=[y_coords[idx]] + list(row))
        
        yscroll = ttk.Scrollbar(tab, orient=tk.VERTICAL, command=tree.yview)
        xscroll = ttk.Scrollbar(tab, orient=tk.HORIZONTAL, command=tree.xview)
        tree.configure(yscroll=yscroll.set, xscroll=xscroll.set)
        
        tree.grid(row=0, column=0, sticky="nsew")
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        
        tab.grid_rowconfigure(0, weight=1)
        tab.grid_columnconfigure(0, weight=1)

    def create_plots_tab(self, notebook, u_data, v_data, uv_data, h, k):
        tab = ttk.Frame(notebook)
        notebook.add(tab, text="Графики")
        
        fig = plt.Figure(figsize=(20, 8))
        
        x = np.linspace(0, 2, u_data.shape[1])
        y = np.linspace(0, 1, u_data.shape[0])
        X, Y = np.meshgrid(x, y, indexing="ij")
        
        ax1 = fig.add_subplot(231, projection='3d')
        ax1.plot_surface(X, Y, u_data.T, cmap='viridis')
        ax1.set_title("Численное решение (u)")
        
        ax2 = fig.add_subplot(232, projection='3d')
        ax2.plot_surface(X, Y, v_data.T, cmap='plasma')
        ax2.set_title("Точное решение (v)")
        
        ax3 = fig.add_subplot(233, projection='3d')
        ax3.plot_surface(X, Y, uv_data.T, cmap='magma')
        ax3.set_title("Разница |u-v|")
        
        ax4 = fig.add_subplot(234)
        contour = ax4.contourf(X, Y, u_data.T, levels=20, cmap='magma')
        fig.colorbar(contour, ax=ax4, shrink=0.7)
        ax4.set_title("u")
        
        ax5 = fig.add_subplot(235)
        contour = ax5.contourf(X, Y, v_data.T, levels=20, cmap='magma')
        fig.colorbar(contour, ax=ax5, shrink=0.7)
        ax5.set_title("v")

        ax6 = fig.add_subplot(236)
        contour = ax6.contourf(X, Y, uv_data.T, levels=20, cmap='magma')
        fig.colorbar(contour, ax=ax6, shrink=0.7)
        ax6.set_title("|u-v|")

        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def create_plots_tab1(self, notebook, u_data, v_data, uv_data, h, k):
        tab = ttk.Frame(notebook)
        notebook.add(tab, text="Графики")
        
        fig = plt.Figure(figsize=(20, 8))
        
        x = np.linspace(0, 2, u_data.shape[1])
        y = np.linspace(0, 1, u_data.shape[0])
        X, Y = np.meshgrid(x, y, indexing="ij")

        x1 = np.linspace(0, 2, v_data.shape[1])
        y1 = np.linspace(0, 1, v_data.shape[0])
        X1, Y1 = np.meshgrid(x1, y1, indexing="ij")
        
        ax1 = fig.add_subplot(231, projection='3d')
        ax1.plot_surface(X, Y, u_data.T, cmap='viridis')
        ax1.set_title("Решение v на сетке (n, m)")
        
        ax2 = fig.add_subplot(232, projection='3d')
        ax2.plot_surface(X1, Y1, v_data.T, cmap='plasma')
        ax2.set_title("Решение v2 на сетке (2n, 2m)")
        
        ax3 = fig.add_subplot(233, projection='3d')
        ax3.plot_surface(X, Y, uv_data.T, cmap='magma')
        ax3.set_title("Разница |v-v2|")
        
        ax4 = fig.add_subplot(234)
        contour = ax4.contourf(X, Y, u_data.T, levels=20, cmap='magma')
        fig.colorbar(contour, ax=ax4, shrink=0.7)
        ax4.set_title("v")
        
        ax5 = fig.add_subplot(235)
        contour = ax5.contourf(X1, Y1, v_data.T, levels=20, cmap='magma')
        fig.colorbar(contour, ax=ax5, shrink=0.7)
        ax5.set_title("v2")

        ax6 = fig.add_subplot(236)
        contour = ax6.contourf(X, Y, uv_data.T, levels=20, cmap='magma')
        fig.colorbar(contour, ax=ax6, shrink=0.7)
        ax6.set_title("|v-v2|")

        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)



if __name__ == "__main__":
    app = MathResultsApp(root)
    root.geometry("1000x800")
    root.mainloop()