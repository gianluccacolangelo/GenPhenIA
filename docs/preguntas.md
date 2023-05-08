Qué tan confiables tienden a ser los fenotipos registrados por los médicos?
Qué lista de genes candidatos son los que mandan a mapear?
Porque ahora entiendo que es más barato secuenciar solo regiones con alta profundidad que mandar a secuenciar todo el genoma.

¿Qué uso como training dataset?

Necesito algún mail para descargar lo de omim. Con HPO me arreglo en cuanto a genotipo-fenotipo. Pero con OMIM terminaría de agregar las enfermedades.

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
		Entonces el modelo básicamente predeciría categorías (genes) en base a si son capitales o específicos... parece muy bruto.

		En cambio... Si realizamos esas métricas pero con respecto al conjunto de fenotipos del paciente, en mi opinión va a ser mejor.
		Esas métricas van a ser nuestros features

