import pickle
import numpy as np

from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

with open('train_featurized_5.p', 'rb') as f:

    train_dataset = pickle.load(f)

train_dataset = np.array(train_dataset)

with open('test_featurized_5.p', 'rb') as f:

    test_dataset = pickle.load(f)

test_dataset = np.array(test_dataset)

print(train_dataset.shape)
print(test_dataset.shape)
#print(dataset[0])

X_train = train_dataset[:,:-1]
y_train = train_dataset[:,-1:]

X_test = test_dataset[:,:-1]
y_test = test_dataset[:,-1:]

#seed = 2018
#test_size = 0.2
#X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)

print('Fitting model')
model = XGBClassifier()
model.fit(X_train, y_train)

print('making predictions')
y_pred = model.predict(X_test)
predictions = [round(value) for value in y_pred]

accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))
