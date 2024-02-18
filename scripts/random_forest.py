##{{{
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np
## }}}


##{{{
diabetes = load_diabetes()
df = pd.DataFrame(data=diabetes.data,columns=diabetes.feature_names)
df['target'] = diabetes.target

X = df.drop('target',axis=1)
y = df['target']

X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=.2,random_state=42)

rf_regressor = RandomForestRegressor(
        n_estimators = 15,
        max_depth = 10,
        min_samples_split = 2,
        max_features=5,
        bootstrap=True,
        random_state=2)

##}}}

##{{{
rf_regressor.fit(X_train,y_train)
##}}}

##{{{
y_pred_train = rf_regressor.predict(X_train)
y_pred_test = rf_regressor.predict(X_test)
mse_test = mean_squared_error(y_test,y_pred_test)
mse_train = mean_squared_error(y_train,y_pred_train)
print(np.sqrt(mse_train))
print(np.sqrt(mse_test))
##}}}


##{{{
importances = rf_regressor.feature_importances_
# Get the standard deviations of the importances
std = np.std([tree.feature_importances_ for tree in rf_regressor.estimators_], axis=0)

# Get feature names
feature_names = X_train.columns

# Create a pandas series with feature importances
feature_importances = pd.Series(importances, index=feature_names).sort_values(ascending=False)

##}}}


##{{{

import matplotlib.pyplot as plt

# Plot the feature importances
plt.figure(figsize=(10, 6))
feature_importances.plot(kind='bar', yerr=std)
plt.title('Feature Importances with Standard Deviation')
plt.xlabel('Features')
plt.ylabel('Importance')
plt.show()
##}}}


##{{{
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score

param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [None, 10, 20, 30],
    'min_samples_split': [2, 4, 6]
}


grid_search = GridSearchCV(estimator=rf_regressor, param_grid=param_grid, cv=5, scoring='neg_mean_squared_error')



grid_search.fit(X_train, y_train)



# Print the best parameters found
print("Best parameters found: ", grid_search.best_params_)

# Evaluate on the test set
best_rf = grid_search.best_estimator_
y_pred = best_rf.predict(X_test)

# Calculate RMSE and R-squared
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
r2 = r2_score(y_test, y_pred)

print(f"RMSE: {rmse}")
print(f"R^2: {r2}")
## }}}


## {{{

simulated_set = pd.read_csv("simulated_dataset.csv")
## }}}
