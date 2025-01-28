pub use model::TFLiteModelData;

/// Pileup model to use
///
/// This reflects the final pileup model selection. It allows us to keep any
/// option synthesis logic out of the model init.
///
#[derive(Default)]
pub enum ModelSelection {
    /// This should only be selected for count mode
    #[default]
    None,

    /// Built-in model
    Builtin,

    /// Override built-in models with a tflite model file, intended for development only
    OverrideFilename(String),
}

#[cfg(not(feature = "tflite"))]
mod model {
    // Create a mock tflite object to compile against when tflite feature is disabled
    use super::*;
    use std::marker::PhantomData;

    pub struct TFLiteModelData<'a> {
        phantom: PhantomData<&'a i32>,
    }

    impl TFLiteModelData<'_> {
        pub fn new(_: &ModelSelection) -> Self {
            Self {
                phantom: PhantomData,
            }
        }
    }
}

#[cfg(feature = "tflite")]
mod model {
    use super::*;
    use tflite::ops::builtin::BuiltinOpResolver;

    pub struct TFLiteModelData<'a> {
        pub interpreter: tflite::Interpreter<'a, BuiltinOpResolver>,
        pub input_index: tflite::TensorIndex,
        pub output_index: tflite::TensorIndex,
    }

    impl TFLiteModelData<'_> {
        pub fn new(model_selection: &ModelSelection) -> Self {
            let model = {
                let builtin_model_path =
                    std::include_bytes!("models/pileup_calling_model.v1.tflite");
                match model_selection {
                    ModelSelection::None => panic!("No methylation model selected"),
                    ModelSelection::Builtin => {
                        tflite::FlatBufferModel::build_from_buffer(builtin_model_path.to_vec())
                    }
                    ModelSelection::OverrideFilename(x) => {
                        tflite::FlatBufferModel::build_from_file(x)
                    }
                }
                .unwrap()
            };
            let resolver = BuiltinOpResolver::default();
            let builder = tflite::InterpreterBuilder::new(model, resolver).unwrap();

            let mut interpreter = builder.build().unwrap();
            interpreter.allocate_tensors().unwrap();
            let inputs = interpreter.inputs().to_vec();
            assert_eq!(inputs.len(), 1);
            let outputs = interpreter.outputs().to_vec();
            assert_eq!(outputs.len(), 1);

            Self {
                interpreter,
                input_index: inputs[0],
                output_index: outputs[0],
            }
        }
    }
}
