cp -r output/fireball_radius python/train/
cd python
cp -r model_parameters train/
python3 process_train_data.py
python3 optimize_grid_regression_models.py

# note: need to run generate_training_data.sh first
