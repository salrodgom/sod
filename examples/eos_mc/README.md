# SOD Boltzmann Monte Carlo Estimator

Este proyecto complementa el código SOD existente aportando una utilidad que obtiene la energía esperada de sustituciones Si→Ge mediante ponderación de Boltzmann. El programa está pensado para compilarse en el mismo entorno que los módulos Fortran del proyecto original, reutilizando rutinas como `init_energy_calc` y `calculate_structure_energy`.

## Requisitos previos

- **Fortran**: compilador compatible con Fortran 2008 (por ejemplo, `gfortran >= 10`).
- **Proyecto SOD anterior**: se asume que el código del proyecto original (especialmente `src/energy_calc.f90` y dependencias) está disponible para enlazarlo.
- Ficheros de entrada generados por SOD en el directorio de trabajo:
  - `INSOD`
  - `SGO`
  - `n00/ENERGIES`, `n01/OUTSOD`, `n01/ENERGIES`, … (al menos hasta el nivel requerido por la expansión de energía)

## Compilación

1. Sitúate en el directorio raíz del proyecto SOD original.
2. Copia el archivo `src/sod_boltzmann_mc.f90` de este workspace al directorio `src/` del proyecto original (o añade la ruta correspondiente al comando de compilación).
3. Añade el nuevo fichero a tu `Makefile` o comando de compilación, por ejemplo:

```bash
# Ejemplo: añadir el ejecutable al Makefile existente
MC_OBJS = sod_boltzmann_mc.o

sod_boltzmann_mc: $(MC_OBJS) $(COMMON_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
```

Asegúrate de que `$(COMMON_OBJS)` contenga `energy_calc.o` y sus dependencias.

> Para aprovechar el paralelismo opcional, compila con soporte OpenMP (`gfortran -fopenmp ...`). Sin ese flag el programa se ejecutará en modo secuencial incluso si se pasa el argumento `omp`.

## Uso

Una vez compilado (`sod_boltzmann_mc`), ejecuta el programa desde el directorio de una simulación (donde residen `INSOD`, `SGO`, `nXX/` y en particular `n01/OUTSOD`). El programa determina el máximo número posible de sustituciones leyendo los sitios representativos listados en `n01/OUTSOD`:

```bash
./sod_boltzmann_mc [temperatura_K] [max_sustituciones] [muestras] [semilla] [sampler] [omp]
```

- `temperatura_K` (opcional): temperatura en Kelvin para el peso de Boltzmann. Valor por defecto `1000 K`.
- `max_sustituciones` (opcional): límite superior del número de sustituciones evaluadas. Por defecto se analizan todos los valores posibles.
- `muestras` (opcional): número de configuraciones Monte Carlo evaluadas cuando el número de combinaciones supera el umbral de enumeración exacta (200000). Por defecto `5000`.
- `semilla` (opcional): semilla entera para el generador aleatorio. Si se omite, se usa el reloj del sistema.
- `sampler` (opcional): modo de muestreo a emplear cuando se activa Monte Carlo. Acepta `uniform` (muestras independientes sin correlación, valor por defecto) o `metropolis` (cadena de Markov con aceptación Metropolis-Hastings que lleva registro de los intentos rechazados).
- `omp` (opcional): habilita o deshabilita el uso de OpenMP en los cálculos estadísticos. Usa `omp`/`1`/`true` para activarlo (requiere compilar con `-fopenmp`) o `noomp`/`0`/`false` para forzar la ejecución secuencial.

El programa escribe en la salida estándar, para cada número de sustituciones, el número de combinaciones analizadas, si se recurrió a muestreo Monte Carlo, la energía mínima y la energía esperada (ponderada por Boltzmann), junto con la desviación estándar y la probabilidad Boltzmann de la configuración mínima.

### Enumeración vs Muestreo

Para mantener tiempos de ejecución razonables, el programa enumera exhaustivamente todas las combinaciones siempre que el número total no supere `200000`. Cuando se rebasa este umbral se cambia automáticamente a un muestreo Monte Carlo. El modo de muestreo depende del argumento `sampler`:

- `uniform`: selecciona configuraciones al azar sin memoria, equivalente a escoger subconjuntos aleatorios independientes.
- `metropolis`: ejecuta una caminata Metropolis-Hastings que propone intercambios singulares y acepta o rechaza configuraciones de acuerdo con el criterio de Boltzmann, contabilizando el número de intentos rechazados.

En ambos modos el número de evaluaciones corresponde al tercer argumento (`muestras`). Si OpenMP está disponible y habilitado, el cálculo de promedios y sumas de pesos se paraleliza automáticamente. En los casos muestreados, la salida advierte explícitamente que no se cubrieron todas las combinaciones.