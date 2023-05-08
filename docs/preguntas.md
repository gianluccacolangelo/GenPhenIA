Qué tan confiables tienden a ser los fenotipos registrados por los médicos?
Qué lista de genes candidatos son los que mandan a mapear?
Porque ahora entiendo que es más barato secuenciar solo regiones con alta profundidad que mandar a secuenciar todo el genoma.

¿Qué uso como training dataset?

Necesito algún mail para descargar lo de omim. Con HPO me arreglo en cuanto a genotipo-fenotipo. Pero con OMIM terminaría de agregar las enfermedades.

Modelo de regresión logística: ¿Qué queremos clasificar? Si son enfermedades, tendríamos un modelo con más de 4mil clases (enfermedades) y métricas muy pobres como la promiscuidad o no de un gen y un fenotipo para determinarlo. Pero si en cambio hacemos para cada enfermedad (con su conjunto de fenotipos) vamos a tener 4mil modelos diferentes, uno para cada enfermedad. La otra también, es para cada paciente, le haces un modelo, que es para solo ese paciente con sus conjuntos de fenotipos, y ahí determinas si un dado gen causa o no causa ese conjunto de fenotipos.
