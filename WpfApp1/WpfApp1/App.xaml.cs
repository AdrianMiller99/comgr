using System.Configuration;
using System.Data;
using System.Windows;
using System.Numerics;
using System.Windows.Media;

namespace WpfApp1;

/// <summary>
/// Interaction logic for App.xaml
/// </summary>
public partial class App : Application
{
    protected override void OnStartup(StartupEventArgs e)
    {
        base.OnStartup(e);

        // Example usage of System.Numerics
        Vector3 vector = new Vector3(1, 2, 3);
        Vector3 anotherVector = new Vector3(4, 5, 6);
        Vector3 result = Vector3.Add(vector, anotherVector);
    }
}