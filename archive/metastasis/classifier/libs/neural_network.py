#!python3
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.utils.np_utils import to_categorical


def neural_network(train, train_labels,
                   hidden_nodes=64, epochs=100):
    """Use Keras Sequential model in the backend.

    Parameters
    ----------
        train: array
            Input data from function `split_dataset`.
        train_labels: array
            Labels from function `split_dataset`.
        hidden_layer: int, optional
            Number of nodes in a hidden layer
        epochs: int, optional.
            Number of one forward pass and one backward pass
            of all the training examples.

    Return
    ------
        A fitted neural network model.
    """
    train_labels = to_categorical(train_labels)
    model = Sequential()
    # number of hidden layer units
    model.add(Dense(hidden_nodes,
                    input_dim=train.shape[1], activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(hidden_nodes, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(4, activation='softmax'))  # Four classes
    model.compile(optimizer='rmsprop',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
    model.summary()
    # TODO: tune hyperparameters
    model.fit(train, train_labels, epochs=epochs, validation_split=0.2)
    return model


def evaluate_nn(model, test, test_labels):
    """Get loss and accuracy of model.

    Parameters
    ----------
        model: keras model
        The neural network model trained
        test: array
            Test set
        test_labels: array
            True test set results
    """
    test_labels = to_categorical(test_labels)
    loss, accuracy = model.evaluate(test, test_labels)
    accuracy *= 100
    print(
        f'Test set result:\nLoss {loss:{1}.{4}} with accuracy {accuracy:{1}.{4}}%.')
    return loss, accuracy / 100
