from deepsort import DeepSortPredictor
# define the model
model = DeepSortPredictor(species='human',
                          tissue='Blood')
# use the trained model to predict
test_files = 'Data/cvid.c2_data.csv'
#['/path/to/human_brain_test_data_1.csv', '/path/to/human_brain_test_data_2.csv']
#for test_file in test_files:
model.predict(test_file, save_path='Output')
