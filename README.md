# Modelo depredador-presa (Lotka-Volterra con capacidad de carga)

Notebook de Jupyter que simula la interacción entre una población de **presas** y una de **depredadores**. Resuelve el sistema de ecuaciones diferenciales de forma numérica, calcula los **puntos de equilibrio** con sus **autovalores** y genera las gráficas de **series de tiempo** y **retrato de fase**.

**Autores:** David Steven Garzón Rodríguez y Carlos Arturo Acevedo Santiago

## Modelo matemático

```
dR/dt = α·R·(1 − R/K) − β·R·F      (presas)
dF/dt = −γ·F + δ·R·F               (depredadores)
```

| Parámetro | Significado | Valor por defecto |
|---|---|---|
| `alpha` (α) | Tasa de crecimiento de las presas | 2.0 |
| `beta` (β) | Tasa de depredación | 1.5 |
| `gamma` (γ) | Mortalidad de los depredadores | 1.0 |
| `delta` (δ) | Eficiencia de conversión presa → depredador | 0.5 |
| `K` | Capacidad de carga del entorno | 2.5 |

## Requisitos

- Python 3.9 o superior (el proyecto fue desarrollado con Python 3.14.7).
- Las dependencias listadas en el archivo `requirements.txt` (se instalan con un solo comando, ver más abajo).

## Instalación

1. **Descarga el proyecto**:

   ```bash
   git clone https://github.com/dsgarzonrodriguez/Sistema-depredador-presa.git
   ```

   O descarga el ZIP desde GitHub y descomprímelo.

2. **Crea un entorno virtual** (recomendado):

   ```bash
   python -m venv .env
   ```

   Actívalo:

   - Windows: `.env\Scripts\activate`
   - Linux / macOS: `source .env/bin/activate`

3. **Instala las dependencias** desde `requirements.txt`:

   ```bash
   pip install -r requirements.txt
   ```

## Uso

1. Abrelo en visual studio code (con la extensión *Jupyter*) o puedes subirlo a google colabs

2. Ejecuta todas las celdas en orden: menú **Run → Run All Cells**, o con `Shift + Enter` celda por celda.

> **Importante:** la imagen del modelo se carga desde `resources/Modelo-depredador-presa.png`. Mantén la carpeta `resources/` junto al notebook, o la imagen no se mostrará.

## Estructura del proyecto

```
.
├── Modelo_depredador_presa-David-Carlos.ipynb   # Código y análisis
├── resources/
│   └── Modelo-depredador-presa.png              # Imagen del modelo
├── requirements.txt                             # Dependencias del proyecto
└── README.md
```

## Cómo usar la clase `Escenario`

Cada escenario es un objeto independiente con sus propios parámetros:

```python
# 1. Crear el escenario con los parámetros del modelo
modelo = Escenario(alpha=2.0, beta=1.5, gamma=1.0, delta=0.8, K=2.5)

# 2. Simular con poblaciones iniciales, tiempo final y paso
modelo.simular(poblacion_presas_inicial=1.0,
               poblacion_depredador_inicial=0.5,
               tiempo_final=16,
               paso=0.01)

# 3. Calcular equilibrios, Jacobiano y autovalores
resultados = modelo.calcular_puntos_equilibrio()

# 4. Graficar
modelo.graficar_series()   # Población vs. tiempo
modelo.retrato_de_fase()   # Presas vs. depredadores (usar dentro de una figura de plt)
```

### Métodos disponibles

| Método | Qué hace |
|---|---|
| `simular(...)` | Resuelve las ecuaciones con RK45 y devuelve `(tiempo, presas, depredadores)`. |
| `calcular_puntos_equilibrio()` | Devuelve, para cada equilibrio, el punto, el Jacobiano y los autovalores. |
| `graficar_series()` | Muestra la evolución de ambas poblaciones en el tiempo. |
| `retrato_de_fase()` | Dibuja la trayectoria en el plano presa-depredador, con flechas de dirección. No crea la figura ni llama a `plt.show()`, así que puedes superponer varios escenarios. |

### Superponer varios escenarios en un retrato de fase

```python
plt.figure(figsize=(10, 5))
caso1.retrato_de_fase()
caso2.retrato_de_fase()
plt.xlabel("Presas (R)")
plt.ylabel("Depredadores (F)")
plt.legend()
plt.grid(True)
plt.show()
```

## Las ecuaciones en el código

El modelo matemático está escrito **dos veces** dentro de la clase `Escenario`, porque cada versión cumple una función distinta.

### 1. Versión simbólica (SymPy), en `__init__`

Se usa para calcular los puntos de equilibrio y el Jacobiano.

```python
# Símbolos algebraicos
self.presas_symb, self.depredador_symb = sp.symbols('R F', real=True)

# Ecuaciones del sistema
self.ecuacion_presas = self.alpha * self.presas_symb * (1 - self.presas_symb / self.K) - self.beta * self.presas_symb * self.depredador_symb
self.ecuacion_depredador = -self.gamma * self.depredador_symb + self.delta * self.presas_symb * self.depredador_symb
```

### 2. Versión numérica (SciPy), en `sistema_edo`

Se usa para simular la evolución en el tiempo con `solve_ivp`.

```python
def sistema_edo(self, tiempo, estado):
    R, F = estado
    dRdt = self.alpha * R * (1 - R / self.K) - self.beta * R * F
    dFdt = -self.gamma * F + self.delta * R * F
    return [dRdt, dFdt]
```

> **Regla de oro:** si cambias una ecuación, cámbiala **en los dos lugares**. Si solo modificas una, la simulación y el análisis de equilibrio describirán modelos distintos.

## Adaptar las ecuaciones a otros modelos

Las ecuaciones son independientes del resto del código: la simulación, el cálculo de equilibrios, el Jacobiano, los autovalores y las gráficas funcionan con cualquier sistema de **dos variables**. Para usar otro modelo solo necesitas editar las dos versiones y, si hace falta, los parámetros del constructor.

### Ejemplo 1: Lotka-Volterra clásico (sin capacidad de carga)

Se elimina el término `(1 − R/K)`:

```python
# En __init__ (simbólica)
self.ecuacion_presas = self.alpha * self.presas_symb - self.beta * self.presas_symb * self.depredador_symb

# En sistema_edo (numérica)
dRdt = self.alpha * R - self.beta * R * F
```

La ecuación de los depredadores no cambia. Este modelo produce ciclos cerrados en lugar de convergir a un punto.

### Ejemplo 2: Respuesta funcional de Holling tipo II

Se agrega un parámetro nuevo `h` (tiempo de manipulación) que satura la caza:

```python
# En el constructor: def __init__(self, ..., K=2.5, h=0.5):
self.h = h

# En __init__ (simbólica)
R, F = self.presas_symb, self.depredador_symb
self.ecuacion_presas = self.alpha * R * (1 - R / self.K) - self.beta * R * F / (1 + self.h * R)
self.ecuacion_depredador = -self.gamma * F + self.delta * R * F / (1 + self.h * R)

# En sistema_edo (numérica)
dRdt = self.alpha * R * (1 - R / self.K) - self.beta * R * F / (1 + self.h * R)
dFdt = -self.gamma * F + self.delta * R * F / (1 + self.h * R)
```

Con ecuaciones más complejas, `sp.solve` puede tardar más o no encontrar todas las soluciones de forma exacta.

### Ejemplo 3: Otros sistemas de dos variables

La misma estructura sirve para otros modelos de dos poblaciones o magnitudes, como competencia entre especies o modelos epidemiológicos de dos variables. Sustituye las ecuaciones, ajusta los nombres de los parámetros y las etiquetas de las gráficas.

### Si tu modelo tiene más de dos variables

Hay que modificar más partes del código, porque hoy está escrito para dos variables:

1. En `sp.symbols(...)`: agrega el nuevo símbolo.
2. En `sistema_edo`: desempaqueta el nuevo valor de `estado`, calcula su derivada y devuélvela en la lista.
3. En `simular`: agrega la condición inicial en `y0` y extrae la nueva fila de `solucion.y`.
4. En `calcular_puntos_equilibrio`: incluye la nueva ecuación y el nuevo símbolo en `sp.solve` y en `.jacobian(...)`.
5. En las gráficas: agrega la nueva serie. El retrato de fase en 2D solo muestra dos variables a la vez.

### Qué más puedes cambiar sin tocar las ecuaciones

| Qué | Dónde |
|---|---|
| Valores de los parámetros | Al crear el escenario: `Escenario(alpha=..., beta=..., ...)` |
| Poblaciones iniciales, duración y paso | Al llamar a `simular(...)` |
| Método de integración | El argumento `method` de `solve_ivp` (por ejemplo `'RK23'`, `'Radau'` o `'LSODA'`) |

## Resultados incluidos

Con α=2, β=1.5, γ=1, δ=0.8 y K=2.5, el notebook encuentra tres puntos de equilibrio:

| Punto (R, F) | Autovalores | Tipo |
|---|---|---|
| (0, 0) | 2 y −1 | Punto silla (inestable) |
| (2.5, 0) | −2 y 1 | Punto silla (inestable) |
| (1.25, 0.67) | −0.5 ± 0.87i | Foco estable (las poblaciones oscilan y convergen) |

El notebook también compara cuatro condiciones iniciales distintas. El `caso4` parte justo del equilibrio coexistente, así que su trayectoria se queda quieta.

## Problemas frecuentes

| Problema | Solución |
|---|---|
| `ModuleNotFoundError: No module named 'sympy'` (u otra librería) | Ejecuta `pip install -r requirements.txt` con el entorno virtual activado. |
| La imagen del modelo no aparece | Verifica que exista `resources/Modelo-depredador-presa.png` junto al notebook. |
| Las gráficas no se muestran | Ejecuta la celda de imports primero y usa un entorno con Jupyter (no un script sin interfaz gráfica). |
| El notebook usa un kernel que no existe en tu equipo | En Jupyter: **Kernel → Change Kernel** y elige tu entorno de Python. |
