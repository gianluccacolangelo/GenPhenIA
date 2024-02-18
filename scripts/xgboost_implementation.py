##{{{
import xgboost as xgb
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
##}}}


##{{{
# Load the diabetes dataset
diabetes_data = load_diabetes()
X, y = diabetes_data.data, diabetes_data.target
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
##}}}


##{{{


xg_reg = xgb.XGBRegressor(objective='reg:squarederror', colsample_bytree=1, learning_rate=0.1,
                          max_depth=100, alpha=1, n_estimators=100)
# Fit the regressor to the training set
xg_reg.fit(X_train, y_train)

# Predict the labels of the test set
y_pred = xg_reg.predict(X_test)
y_train_pred = xg_reg.predict(X_train)

# Compute and print the mean squared error of the predictions
rmse = mean_squared_error(y_test, y_pred, squared=False)
bias = mean_squared_error(y_train,y_train_pred,squared=False)
print(f"RMSE: {rmse}")
print(f"bias: {bias}")

##}}}
