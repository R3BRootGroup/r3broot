# Class interfaces

## `R3BIOConnector`

The class `R3BIOConnector` is used to declare the input and output data structure in the tasks based on `FairTask`. It's highly recommended to use `R3BIOConnector` instead of `TClonesArray`.

Usage:

```c++
class MyTask : public FairTask
{
  public:
    MyTask() = default;
    auto Init() -> InitStatus override
    {
        points_.init();
        hits_.init();
    }

    void Exec(Option_t* /*option*/) override {
        const auto& points = points_.get();
        // do something about points
        // ...

        // fill the hits
        auto& hits = hits_.get();
        auto& hit = hits.emplace_oback();
    }

  private:
    R3B::InputVectorConnector<R3BNeulandPoint> points_{ "NeulandPoints" };
    R3B::OutputVectorConnector<R3BNeulandPoint> hits_{ "NeulandHits" };
}
```

Additionally there are three built-in types to support data IO with STL containers:

- `InputVectorConnector` and `OutputVectorConnector`: data IO with `std::vector`
- `InputMapConnector` and `OutputMapConnector`: data IO with `std::map`
- `InputHashConnector` and `OutputHashConnector`: data IO with `std::unordered_map`

The `get()` API provides the access to the internal STL data container from `R3BIOConnector`. For `R3BInputConnector`, the return value is a const reference `const T&`. And for `R3BOutputConnector`, the return value is a simple reference `T&`.

Additionally, the `R3BInputConnector` provides iterators to loop through every element in the internal STL container:
```c++
// iterate through points using range-based for loop
for(const auto& point : points_)
{
    std::cout << point.energy << std::endl;
}
```

**Attention:**

To successfully read and write STL container from/to the root file, the corresponding types must be also declared in the `LinkDef.h` file. For example:

```c++
#pragma link C++ class R3BNeulandPoint+;
#pragma link C++ class vector<R3BNeulandPoint>+;
```


## `R3BValueError`

`R3BValueError` is a data structure which contains a value and an error to represent any measurement value with an uncertainty. 

Usage:
```c++
const auto distance = R3B::ValueError{1., 0.2}; // represents 1 +/- 0.2
```

Some basic arithmetic operations are implemented, such as `+`, `*`, `/`, `-`, `-=` and `+=`. For example:

```c++
const auto v1 = R3B::ValueError{1., 0.2};
const auto v2 = R3B::ValueError{3., 0.2};

const auto v3 = v1 * v2;
const auto v4 = v1 / v2;
const auto v5 = v1 + v2;
const auto v6 = v1 - v2;
```
Error values from the above mentioned operations are calculated using Gaussian error propagation. Therefore, it could be expensive in terms of computations.

## `R3BFileSource2`

`R3BFileSource2` is an improved version of `R3BFileSource` with a much cleaner design using the same interfaces. It also supports reading root files containing just a single root tree.

Usage:
```c++
// To add multiple files
auto file_source = std::make_unique<R3BFileSource2>();

for(auto filename : filenames)
{
    file_source->AddFile(std::move(filename));
}
run->SetSource(file_source.release()); // if FairRoot version < 19.0.1
run->SetSource(std::move(file_source)); // if FairRoot version >= 19.0.1
```
To add the source files which only contain root trees:
```c++
auto file_source = std::make_unique<R3BFileSource2>();

file_source->SetInitRunID(999); // set the run iD

file_source->AddFile("filename.root", true);
```
Please make sure that the run ID set to file_source is consistent with the run ID from the parameter file.

### Event processing rate printing

`R3BFileSource2` also has a Event processing rate printer, to the set refresh rate of the printer, use

```c++
file_source->SetEventPrintRefreshRate(2); // 2 Hz
```

## `R3BDataMonitor`

`R3BDataMonitor` is a histogram/graph manager for tasks based on `FairTask`.

Usage:
```c++
class MyTask : public FairTask
{
  public:
    MyTask() = default;
    auto Init() -> InitStatus override
    {
        hist_1d = data_monitor_.add_hist<TH1D>(
            "module_num", "Counts with module ids", total_bar_num, -0.5, total_bar_num + 0.5);

        hist_2d = data_monitor_.add_hist<TH2D>("TimeLVsBar",
                                               "Time_vs_bar_left",
                                               total_bar_num,
                                               0.5,
                                               0.5 + total_bar_num,
                                               BINSIZE_TIME,
                                               HIST_MIN_TIME,
                                               HIST_MAX_TIME);
    }

    void Exec(Option_t* /*option*/) override
    {
        hist_1d->Fill(val1);
        hist_2d->Fill(val2, val3);
    }

    void FinishTask() override
    {
        // Save to the sink file in "DataMonitor/foler_name" folder
        histograms_.save_to_sink(folder_name);
    }

  private:
    R3BDataMonitor data_monitor_;
    TH1D* hist_1d = nullptr;
    TH2D* hist_2d = nullptr;
}
```
With `R3BDataMonitor`, all histograms/graphs will be automatically saved and written to the file in the end of the task.

### Output file location

There are two APIs from `R3BDataMonitor` to save histograms/graphs to a local root file:

- `save_to_sink(foldername)`
  saves all the histograms and graphs to the sink file (obtained from `FairRun::Instance()->GetSink()`). The folder where histograms are located is "DataMonitor/foldername". If the `foldername` is empty or not provided, the location would just be "DataMonitor" folder.

- `save_to_file(filename)`
  saves all histograms and graphs to an independent root file with the file name `filename`. If the `filename` is empty or not provided, the file would be named with the current data and time.

### Canvases for online monitoring

Histograms and graphs can also be grouped into different canvases (ROOT `TCavnas`). This could happen quite often for the online monitoring using the ROOT `THttp` server.

To create a new canvas and add histograms to it:

```c++
auto Init() -> InitStatus override
{
        auto& canvas = data_monitor_.creat_canvas("canvas name", "canvas title", x, y, a, b);
        canva.divide(1, 2);
        hHitEvsBarCosmics_ = canvas.add<1, TH2D>("hHitEvsBarCosmics",
                                                 "HitLevel: Energy vs Bars cosmics",
                                                 bar_numbers,
                                                 -0.5,
                                                 bar_numbers + 0.5,
                                                 ENERGY_BIN_SIZE,
                                                 0,
                                                 ENERGY_MAX);

        hTdiffvsBarCosmics_ = canvas.add<2, TH2D>("hTdiffvsBarCosmics",
                                                  "Tdiff vs Bars cosmics",
                                                  bar_numbers,
                                                  -0.5,
                                                  bar_numbers + 0.5,
                                                  TIME_BIN_SIZE,
                                                  -ENERGY_MAX,
                                                  ENERGY_MAX);
}
```
Again, all the histograms and graphs inside each canvas will be automatically saved into the specified folder.

### Canvas interaction with Root `THttp` server class

If the online analysis is based on Root `THttp` server class, `DataMonitor` provides an additional API to register all owning canvases:

```c++
auto* run = FairRun::Instance();
data_monitor_.register_canvases(run);
```
