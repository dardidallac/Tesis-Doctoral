# ============================================
# Comparación de índices de uso de codones
# (VERSIÓN OPTIMIZADA)
# ============================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import combinations
from scipy.stats import mannwhitneyu

# ============================
# Parámetros
# ============================

CSV_FILE = "/home/mario/Escritorio/boxplot indices/resultados_indices_medicago_final.csv"
COL_WRAP = 3
STAR_FONTSIZE = 10

# ============================
# Funciones auxiliares
# ============================

def p_to_star(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return None


def draw_significance(ax, x1, x2, y, h, star):
    ax.plot([x1, x1, x2, x2],
            [y, y + h, y + h, y],
            lw=0.5, c="black")
    ax.text((x1 + x2) / 2,
            y + h * -0.2,
            star,
            ha="center",
            va="bottom",
            fontsize=STAR_FONTSIZE)


# ============================
# Leer y reestructurar datos (Versión Robusta)
# ============================

# 1. Cargamos el CSV sin ninguna interpretación automática de decimales/miles
df = pd.read_csv(CSV_FILE)

# 2. Reestructuramos a formato largo
df_long = df.melt(
    id_vars=["cds", "Sample"],
    var_name="Index",
    value_name="value"
)

# 3. Función de limpieza total
def fix_numeric_values(x):
    if pd.isna(x): 
        return np.nan
    
    # Convertimos a string y quitamos espacios
    val = str(x).strip()
    
    # CASO CRÍTICO: Si el valor tiene una coma, asumimos que es el decimal.
    # Pero antes, quitamos cualquier punto que actúe como separador de miles.
    # Ejemplo: "1.200,50" -> "1200.50"
    if "," in val:
        val = val.replace(".", "") # Quita punto de miles
        val = val.replace(",", ".") # Cambia coma por punto decimal
    
    try:
        return float(val)
    except ValueError:
        return np.nan

# Aplicamos la limpieza
df_long["value"] = df_long["value"].apply(fix_numeric_values)

# 4. Control de calidad: Filtramos outliers imposibles para GC3 (ej. valores > 100)
# Si tus datos de GC3 están entre 0 y 100, esto evitará que los asteriscos vuelen.
df_long = df_long[df_long["value"] <= 100]

# Eliminamos nulos
df_long = df_long.dropna(subset=["value"])

# VERIFICACIÓN: Imprime los primeros valores para estar seguros
print("--- Verificación de datos procesados ---")
print(df_long.groupby("Index")["value"].describe())

# Ahora sí, convertimos a número real (float)
df_long["value"] = pd.to_numeric(df_long["value"], errors='coerce')

# Eliminamos las filas que no pudieron convertirse (encabezados repetidos o nulos)
df_long = df_long.dropna(subset=["value"])

# ---------------------------------------------------------

samples = sorted(df_long["Sample"].unique())
indices = df_long["Index"].unique()
comparisons = list(combinations(samples, 2))

# ============================
# PRECOMPUTAR datos por grupo
# ============================

grouped = {
    (idx, s): g["value"].values
    for (idx, s), g in df_long.groupby(["Index", "Sample"])
}

# ============================
# PRECOMPUTAR TODOS los tests
# ============================

pvals = {}

for idx in indices:
    for s1, s2 in comparisons:
        v1 = grouped[(idx, s1)]
        v2 = grouped[(idx, s2)]

        stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
        pvals[(idx, s1, s2)] = p

# ============================
# Gráficos
# ============================

sns.set_style("whitegrid")
sns.set_context(
    "paper",
    rc={
        "axes.titlesize": 15,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    }
)

g = sns.catplot(
    data=df_long,
    x="Sample",
    y="value",
    col="Index",
    kind="box",
    col_wrap=COL_WRAP,
    sharey=False,
    height=2.6,
    aspect=0.95,
    linewidth=0.9,
    fliersize=2
)

# ============================
# Anotaciones (solo dibujar)
# ============================

for ax, idx in zip(g.axes.flat, indices):

    sub = df_long[df_long["Index"] == idx]

    y_min = sub["value"].min()
    y_max = sub["value"].max()
    y_range = y_max - y_min

    y_step = y_range * 0.06
    y = y_max + y_step

    for s1, s2 in comparisons:
        p = pvals[(idx, s1, s2)]
        star = p_to_star(p)

        if star:
            x1 = list(samples).index(s1)
            x2 = list(samples).index(s2)

            draw_significance(
                ax,
                x1,
                x2,
                y,
                y_step * 0.3,
                star
            )
            y += y_step

    ax.set_title(idx, fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")

# ============================
# Layout final
# ============================

g.fig.subplots_adjust(
    top=0.88,
    hspace=0.45,
    wspace=0.25
)

g.fig.suptitle(
    "Comparación de índices de uso de codones",
    fontsize=12
)

plt.show()

