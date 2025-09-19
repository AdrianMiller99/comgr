/*// See https://aka.ms/new-console-template for more information

Console.WriteLine("Hello, World!");*/

using System;
using System.Windows.Forms;
using Gradient;

internal static class Program
{
    [STAThread]
    static void Main()
    {
        ApplicationConfiguration.Initialize();
        Application.Run(new GradientForm());
    }
}
