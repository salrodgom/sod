# SOD Ensemble Monte Carlo Estimator

Este proyecto complementa el código SOD existente aportando una utilidad que obtiene la energía esperada de sustituciones Si→Ge mediante ponderación de Boltzmann. El programa está pensado para compilarse en el mismo entorno que los módulos Fortran del proyecto original, reutilizando rutinas como `init_energy_calc` y `calculate_structure_energy`.

> Nota: a partir de esta refactorización el ejecutable principal es `sod_ensemble`, que acepta el modo `mc` para invocar internamente el mismo flujo lógico descrito aquí. El binario `sod_ensemble_mc` sigue generándose como un envoltorio compatibilizador y todos los argumentos se mantienen.

## Requisitos previos

- **Fortran**: compilador compatible con Fortran 2008 (por ejemplo, `gfortran >= 10`).
- **Proyecto SOD anterior**: se asume que el código del proyecto original (especialmente `src/sod_ensemble_energy_calculations.f90` y dependencias) está disponible para enlazarlo.
- Ficheros de entrada generados por SOD en el directorio de trabajo:
  - `INSOD`
  - `SGO`
  - `n00/ENERGIES`, `n01/OUTSOD`, `n01/ENERGIES`, … (al menos hasta el nivel requerido por la expansión de energía)

## Compilación

1. Sitúate en el directorio raíz del proyecto SOD original.
2. Copia el archivo `src/sod_ensemble_mc.f90` de este workspace al directorio `src/` del proyecto original (o añade la ruta correspondiente al comando de compilación).
3. Añade el nuevo fichero a tu `Makefile` o comando de compilación, por ejemplo:

```bash
# Ejemplo: añadir el ejecutable al Makefile existente
MC_OBJS = sod_ensemble_mc.o

sod_ensemble_mc: $(MC_OBJS) $(COMMON_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
```

Asegúrate de que `$(COMMON_OBJS)` contenga `sod_ensemble_energy_calculations.o` y sus dependencias.

> Para aprovechar el paralelismo opcional, compila con soporte OpenMP (`gfortran -fopenmp ...`). Sin ese flag el programa se ejecutará en modo secuencial incluso si se pasa el argumento `omp`.

## Uso

Una vez compilado (`sod_ensemble_mc`), ejecuta el programa desde el directorio de una simulación (donde residen `INSOD`, `SGO`, `nXX/` y en particular `n01/OUTSOD`). El programa determina el máximo número posible de sustituciones leyendo los sitios representativos listados en `n01/OUTSOD`:

```bash
./sod_ensemble_mc -T 1000 -M 12 -C 5000 -s 1234 -a metropolis --omp -N 3:8
```

Parámetros principales (puedes combinar abreviaturas cortas y largas; la sintaxis posicional clásica sigue funcionando para compatibilidad):

- `-T`, `--temperature <K>`: temperatura en Kelvin para los pesos de Boltzmann. Por defecto `1000`.
- `-M`, `--max-substitutions <N>`: máximo de sustituciones a evaluar si no se usa `-N`. Por defecto `-1` (todos los niveles permitidos por `OUTSOD`).
- `-C`, `--samples <N>`: número de configuraciones Monte Carlo evaluadas cuando `C(N,npos)` supera el umbral de enumeración exacta (200000). Por defecto `5000`.
- `-s`, `--seed <valor>`: semilla entera para el generador aleatorio. Si se omite o vale `-1`, se inicializa con `system_clock`.
- `-a`, `--sampler <modo>`: modo de muestreo (`uniform` por defecto o `metropolis`).
- `--omp` / `--no-omp`: fuerza el uso de OpenMP (requiere compilar con `-fopenmp`) o ejecuta siempre de forma secuencial.
- `-N <spec>`: rango o lista de niveles (igual que antes: `-N 5`, `-N 3:8`, `-N 12,30,45`).
- `--osda-gin`, `--no-osda-gin`, `--force-mc`, `--parallel-lists`: opciones adicionales idénticas a la versión anterior.

Ejemplo con la sintaxis heredada (todos los argumentos posicionales):

```bash
./sod_ensemble_mc 1000 12 5000 1234 metropolis omp -N 3:8
```
El programa escribe en la salida estándar, para cada número de sustituciones, el número de combinaciones analizadas, si se recurrió a muestreo Monte Carlo, la energía mínima y la energía esperada (ponderada por Boltzmann), junto con la desviación estándar y la probabilidad Boltzmann de la configuración mínima.

### Enumeración vs Muestreo

Para mantener tiempos de ejecución razonables, el programa enumera exhaustivamente todas las combinaciones siempre que el número total no supere `200000`. Cuando se rebasa este umbral se cambia automáticamente a un muestreo Monte Carlo. El modo de muestreo depende del argumento `sampler`:

- `uniform`: selecciona configuraciones al azar sin memoria, equivalente a escoger subconjuntos aleatorios independientes.
- `metropolis`: ejecuta una caminata Metropolis-Hastings que propone intercambios singulares y acepta o rechaza configuraciones de acuerdo con el criterio de Boltzmann, contabilizando el número de intentos rechazados.

En ambos modos el número de evaluaciones corresponde al tercer argumento (`muestras`). Si OpenMP está disponible y habilitado, el cálculo de promedios y sumas de pesos se paraleliza automáticamente. En los casos muestreados, la salida advierte explícitamente que no se cubrieron todas las combinaciones.