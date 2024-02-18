Qué tan confiables tienden a ser los fenotipos registrados por los médicos?
Qué lista de genes candidatos son los que mandan a mapear?
Porque ahora entiendo que es más barato secuenciar solo regiones con alta profundidad que mandar a secuenciar todo el genoma.

¿Qué uso como training dataset?

Necesito algún mail para descargar lo de omim. Con HPO me arreglo en cuanto a genotipo-fenotipo. Pero con OMIM terminaría de agregar las enfermedades.

En omim: las enfermedades que están asociadas a más de un gen, significa que ocurren cuando ocurren mutaciones en los 3 genes o que cualquiera de los 3 independientemente causa la enfermedad? Pregunta filosófica, no sería una enfermedad diferente?

## Regresión logística:
Primero tenemos que pensar qué categoría nos interesa predecir:

* Un gen específico entre 4000 (modelo logístico general para cualquier conj fenotipos, absurdo de computar)
* Un conjunto de genes (pre-clusterizado) entre ~70. (modelo logístico general razonable de computar, pero dudo de las métricas generales)
* Una lista de genes entre 200 (filtrando de los 4000 según los fenotipos que son dados) (modelo logístico específico del paciente o conj de fenotipos. Mejor rendimiento en mi opinión, pero lo negativo es que es solamente para un conj de fenotipos)



Mis variables (features) son los fenotipos
Mis clases a asignar son los genes.

Problema: Hay más de 4mil genes, es inviable prácticamente hacer una regresión logística multiclase de 4000 clases.
Posible solución: Clusterizar genes.
 	- ¿Cómo?
		K-means
		PCA
		gen2vec
	- Otra solución: trabajar con una lista de genes candidatos y poner métricas (que van a servir de features) particulares para cada paciente.


¿Pero cómo meto los fenotipos en el modelo?
 	- Aquí viene el problema de las métricas:
		Si hacemos unas métricas generales de toda la db (como la relación fenotipo:ngenes y gen:nfenotipos).
		Entonces el modelo básicamente predeciría categorías (genes) en base a si son capitales o específicos... parece muy bruto, además de que todos los genes que estén linkeados a una enfermedad van a tender a tener estos dos parámetros altos... Entonces no podrías diferenciarlos con esas métricas.

		En cambio... Si realizamos esas métricas pero con respecto al conjunto de fenotipos del paciente, en mi opinión va a ser mejor.
		Esas métricas van a ser nuestros features



Tercera pasada en limpio:

Una tercera opción, sería ir tomando enf por enf.

Enfermedad 1. Tomamos todos los fenotipos registrados a la misma, y luego tomamos la unión de los genes que causan esos fenotipos y calculamos, la métrica para el gen que lo causa.
Luego a la enf 2 y así. Entonces tenemos una tabla que tiene

capitalidad | especificidad | enferemedad


viene un paciente y para cada uno de los genes posibles, le calculamos la capitalidad y especificidad, con lo cual en la regresión logística vemos si cae adentro o no.

Pero esto tiene el problema de que la métricas que tenemos miden capitalidad y especificidad, y en el entrenamiento, todos los genes asociados a una enfermedad van a tender a tener estos parámetros altos, entonces no se van a diferenciar... Lo que sí funcionaría sería de una enfermedad, agarrar todos los fenotipos y calcular la métrica para todos los genes posibles. Ahí podríamos ver diferencias de los genes que no lo causan vs lo que sí lo causan.

Está bueno porque para una misma enfermedad quizás podemos tener genes diferentes, y los parámetros de capitalidad y especificidad de esos genes va a cambiar.


Un approach diferente desde el principio:

Pasar los fenotipos de cada enfermedad a un vector de 0 y 1. Donde cada posición es un fenotipo único y si tiene un 1 es porque tiene ese fenotipo. Entonces cuando viene un paciente con fenotipos, se convierte eso a una matriz y se mide la distancia coseno con respecto a todas las enfermedades y nos devuelve las x más cortas. También podríamos usar una regresión logística pero no sé si es lo más eficiente para este approach generar 4000 features con otras 4000 clases diferentes para predecir.

El problema del enfoque de considerar el "vector" de fenotipos es que hay una sola muestra de fenotipos por enfermedad o por gen. De modo que dificulta el entrenamiento y puede caer en el overfitting. Además, un paciente cuando le identifican fenotipos, no tiene necesariamente todos los que están asociados a la enfermedad incógnita en la base de datos. Así, necesitamos un set de entrenamiento. E
